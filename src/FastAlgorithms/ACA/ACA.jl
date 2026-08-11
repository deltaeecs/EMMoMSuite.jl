module ACA

using LinearAlgebra

export LowRankBlock, aca, recompress!, compression_stats

"""
    LowRankBlock{T}

低秩块：`Z ≈ U * transpose(V)`（转置约定，无共轭），其中 `U` 为 m×k、`V` 为 n×k。
对复数对称阻抗矩阵（EFIE 等，`Z = transpose(Z)`），该约定使转置块应用无需共轭。
"""
struct LowRankBlock{T}
    U::Matrix{T}
    V::Matrix{T}
end

Base.size(B::LowRankBlock) = (size(B.U, 1), size(B.V, 1))
Base.size(B::LowRankBlock, i::Int) = size(B)[i]

"""
    aca(Z::AbstractMatrix; tol=1e-4, maxrank=min(size(Z)...), recompress=true)

对稠密块 `Z` 做部分主元自适应交叉近似（Gibson Ch9 Algorithm 6，复数采用转置约定），
返回 `LowRankBlock(U, V)` 满足 `Z ≈ U * transpose(V)`。

# 参数
- `tol`：ACA 收敛容差（相对 Frobenius 范数估计），典型值 `1e-4`~`1e-6`。
- `maxrank`：最大秩上限（默认 `min(m,n)`）。
- `recompress`：是否在迭代后执行 QR/SVD 再压缩（`τ_SVD = 10*tol`）。
"""
function aca(
    Z::AbstractMatrix{T};
    tol::Real = 1e-4,
    maxrank::Int = min(size(Z)...),
    recompress::Bool = true,
) where {T}
    m, n = size(Z)
    getrow(i) = Z[i, :]
    getcol(j) = Z[:, j]
    return aca(T, getrow, getcol, m, n; tol = tol, maxrank = maxrank, recompress = recompress)
end

"""
    aca(::Type{T}, getrow, getcol, m, n; tol=1e-4, maxrank=min(m,n), recompress=true)

自适应交叉近似：通过按行/列采样构造 `Z ≈ U * transpose(V)`。收敛判据为
`‖u_k‖‖v_k‖ ≤ tol * ‖Z̃‖_F`，其中 `‖Z̃‖_F²` 用递推估计（Algorithm 6 step 12）。

行残差 `R(I_k,:) = Z(I_k,:) − Σ_l (u_l)_{I_k} v_l`，列残差
`R(:,J_k) = Z(:,J_k) − Σ_l (v_l)_{J_k} u_l`（均无共轭，转置约定）。
"""
function aca(
    ::Type{T},
    getrow::Function,
    getcol::Function,
    m::Int,
    n::Int;
    tol::Real = 1e-4,
    maxrank::Int = min(m, n),
    recompress::Bool = true,
) where {T}
    # 1. 寻找首个非零行（Gibson 9.2.1.1：避免零行/零块导致失败）
    Iprev = 0
    for i in 1:m
        if norm(getrow(i)) > 0
            Iprev = i
            break
        end
    end
    Iprev == 0 && return LowRankBlock(zeros(T, m, 0), zeros(T, n, 0))

    U = Matrix{T}(undef, m, 0)
    V = Matrix{T}(undef, n, 0)
    rows_used = falses(m)
    cols_used = falses(n)
    normZ2 = 0.0
    eps_t = eps(float(real(T)))
    k = 0

    while k < maxrank
        # 6. 更新第 Iprev 行残差：R(I,:) = Z(I,:) - Σ_l u_l[I] * v_l
        Rrow = k == 0 ? getrow(Iprev) : getrow(Iprev) .- (V * U[Iprev, :])

        # 7. 早期终止：剩余列残差近似为 0
        Jk = 0
        best = 0.0
        for j in 1:n
            cols_used[j] && continue
            a = abs(Rrow[j])
            if a > best
                best = a
                Jk = j
            end
        end
        best < eps_t * max(norm(Rrow), one(real(T))) && break

        # 8-9. 主元列 Jk 与 v_k = R(I_k,:) / R(I_k,J_k)
        pivot = Rrow[Jk]
        v = Rrow ./ pivot

        # 10-11. 列残差 R(:,J_k) = Z(:,J_k) - Σ_l v_l[J_k] * u_l，u_k = R(:,J_k)
        Rcol = k == 0 ? getcol(Jk) : getcol(Jk) .- (U * V[Jk, :])
        u = Rcol

        # 12. Frobenius 范数递推估计：
        #     ‖Z̃(k)‖² = ‖Z̃(k-1)‖² + 2 Re(Σ_j (u_j^H u_k)(v_j^H v_k)) + ‖u_k‖²‖v_k‖²
        s = 0.0
        for l in 1:k
            s += 2 * real(dot(U[:, l], u) * dot(V[:, l], v))
        end
        normZ2 += s + real(dot(u, u) * dot(v, v))

        U = hcat(U, u)
        V = hcat(V, v)
        rows_used[Iprev] = true
        cols_used[Jk] = true
        k += 1

        # 13. 收敛判据
        if norm(u) * norm(v) <= tol * sqrt(max(normZ2, 0.0))
            break
        end

        # 14. 下一主元行：|R(i, J_k)| 最大且未使用
        Iprev = 0
        best = 0.0
        for i in 1:m
            rows_used[i] && continue
            a = abs(Rcol[i])
            if a > best
                best = a
                Iprev = i
            end
        end
        Iprev == 0 && break
    end

    B = LowRankBlock(U, V)
    recompress && (B = recompress!(B; tol = 10 * tol))
    return B
end

"""
    recompress!(B::LowRankBlock; tol=1e-3)

QR/SVD 再压缩（Gibson 9.2.2，转置约定）：`Z = Q_u R_u R_vᵀ Q_vᵀ`，对
`R_u * R_vᵀ` 做 SVD，按 `σ ≥ tol * σ_max` 截断。通常可再获得 20-30% 额外压缩，
`τ_SVD ≈ 10 * τ_ACA`。
"""
function recompress!(B::LowRankBlock{T}; tol::Real = 1e-3) where {T}
    U, V = B.U, B.V
    k = size(U, 2)
    k == 0 && return B

    Fu = qr(U)
    Fv = qr(V)
    Ru = Matrix(Fu.R)
    Rv = Matrix(Fv.R)
    S = Ru * transpose(Rv)
    F = svd(S)
    σmax = F.S[1]
    r = count(s -> s >= tol * σmax, F.S)
    r < 1 && (r = 1)

    U_new = Matrix(Fu.Q) * (F.U[:, 1:r] * Diagonal(F.S[1:r]))
    V_new = Matrix(Fv.Q) * conj(F.V[:, 1:r])
    return LowRankBlock(U_new, V_new)
end

"""
    compression_stats(B::LowRankBlock, m, n) -> (rank, ratio)

压缩统计：`ratio = 1 - k*(m+n)/(m*n)`，`k` 为压缩秩。
"""
function compression_stats(B::LowRankBlock, m::Int, n::Int)
    k = size(B.U, 2)
    return k, 1 - k * (m + n) / (m * n)
end

end # module ACA
