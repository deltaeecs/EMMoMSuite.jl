module CFIEModule

using ..CoreModule
using ..Geometry
using ..BasisFunctions
using ..EFIEModule
using ..MFIEModule
using ..Impedance
using LinearAlgebra
using StaticArrays
using Base.Threads

import ..CoreModule: assemble_impedance_matrix
import ..EFIEModule: efie_interaction!
import ..MFIEModule: mfie_interaction!

export CFIE, assemble_impedance_matrix

"""
    CFIE{FT, CT} <: AbstractIntegralOperator

Combined Field Integral Equation (CFIE) operator.

A linear combination of EFIE and MFIE used to eliminate internal resonance problems that plague EFIE and MFIE for closed structures at specific frequencies.

# Mathematical Formulation

The CFIE is defined as:
```math
\\text{CFIE} = \\alpha \\cdot \\text{EFIE} + (1-\\alpha) \\eta \\cdot \\text{MFIE}
```

The impedance matrix is constructed as:
```math
\\mathbf{Z}_{CFIE} = \\alpha \\mathbf{Z}_{EFIE} + (1-\\alpha) \\eta \\mathbf{Z}_{MFIE}
```

where:
- \$\\alpha\$ is the weighting factor (typically 0.5).
- \$\\eta\$ is the intrinsic impedance of the medium (used to balance units).

# Fields
- `freq`: Operating frequency.
- `alpha`: Weighting factor \$\\alpha\$ (0 to 1).
- `efie`: Underlying `EFIE` operator.
- `mfie`: Underlying `MFIE` operator.
"""
struct CFIE{FT<:AbstractFloat,CT<:Complex} <: AbstractIntegralOperator
    freq::FT
    alpha::FT
    efie::EFIE{FT,CT}
    mfie::MFIE{FT,CT}
end

function CFIE(freq::FT, alpha::FT = 0.5) where {FT}
    efie = EFIE(freq)
    mfie = MFIE(freq)
    return CFIE{FT,Complex{FT}}(freq, alpha, efie, mfie)
end

"""
    assemble_impedance_matrix(cfie::CFIE, basis::RWGBasis)

仿照 Legacy `impedancemat4CFIE4PEC`：单次循环、直接填充唯一的 Z_cfie 矩阵。

- 只分配 **一个** N×N 矩阵（Legacy 风格，原先需要 Z_efie + Z_mfie 两个）
- EFIE 与 MFIE 的远场均使用同一套 4 点 GQ，quadrature points 预计算一次共享
- 每个线程使用栈上 3×3 局部缓冲（Z_e, Z_m），无内部堆分配
- 在写入全局矩阵前就地合并：Z_cfie = α·Z_efie + (1-α)·Z_mfie

注：Z_mfie 已在 `mfie_interaction!` 内乘以 η，此处系数 (1-α) 直接使用。
"""
function assemble_impedance_matrix(cfie::CFIE{FT,CT}, basis::RWGBasis{IT,FT}) where {IT,FT,CT}
    N = num_basis(basis)
    Z = zeros(CT, N, N)           # 唯一的 N×N 分配

    mesh = basis.mesh
    nt = num_elements(mesh)
    efie = cfie.efie
    mfie = cfie.mfie
    α = cfie.alpha
    β = FT(1) - α

    # ── 步骤 1：预计算所有三角形信息 ────────────────────────────────────
    tris_info = Vector{TriangleInfo{IT,FT}}(undef, nt)
    Threads.@threads for t = 1:nt
        tris_info[t] = get_triangle_info(mesh, basis, t)
    end

    # ── 步骤 2：预计算远场 quadrature points（EFIE 与 MFIE 共享同一套 4 点 GQ）
    # EFIE.gq_far 和 MFIE.gq_info 均为 GaussQuadratureInfo(:Triangle, 4, FT)，坐标相同。
    gq = efie.gq_far
    N_gq = length(gq.weight)
    quad_points = Vector{SVector{N_gq,SVector{3,FT}}}(undef, nt)
    Threads.@threads for t = 1:nt
        v = mesh.triangles[:, t]
        v1 = SVector{3,FT}(mesh.node[:, v[1]])
        v2 = SVector{3,FT}(mesh.node[:, v[2]])
        v3 = SVector{3,FT}(mesh.node[:, v[3]])
        quad_points[t] = SVector{N_gq,SVector{3,FT}}(
            v1 * gq.coordinate[1, i] + v2 * gq.coordinate[2, i] + v3 * gq.coordinate[3, i] for
            i = 1:N_gq
        )
    end

    # ── 步骤 3：并行单循环，直填 Z_cfie ─────────────────────────────────
    n_threads = Threads.nthreads()
    row_locks = [SpinLock() for _ = 1:N]

    Threads.@threads for tid = 1:n_threads
        # 线程私有 3×3 缓冲（栈分配，内层循环零堆开销）
        Z_e = zeros(CT, 3, 3)   # EFIE 局部结果
        Z_m = zeros(CT, 3, 3)   # MFIE 局部结果

        for t_test = tid:n_threads:nt
            tri_test = tris_info[t_test]

            for t_src = 1:nt          # MFIE 不对称，遍历全部源三角形
                tri_src = tris_info[t_src]

                fill!(Z_e, zero(CT))
                fill!(Z_m, zero(CT))

                # EFIE 贡献（内部已处理奇异/近奇异，并乘以 efie.factor）
                efie_interaction!(Z_e, efie, tri_test, tri_src, quad_points)
                # MFIE 贡献（内部已乘以 η）
                mfie_interaction!(Z_m, mfie, tri_test, tri_src, quad_points)

                # 就地合并：Z_cfie_local = α·Z_efie + (1-α)·Z_mfie
                @inbounds for j = 1:3, i = 1:3
                    Z_e[i, j] = α * Z_e[i, j] + β * Z_m[i, j]
                end

                # 写入全局矩阵（逐行加锁）
                @inbounds for i = 1:3
                    row_idx = tri_test.inBfsID[i]
                    row_idx == 0 && continue
                    sign_t = tri_test.bfsSign[i]

                    lock(row_locks[row_idx])
                    @inbounds for j = 1:3
                        col_idx = tri_src.inBfsID[j]
                        if col_idx != 0
                            Z[row_idx, col_idx] += Z_e[i, j] * sign_t * tri_src.bfsSign[j]
                        end
                    end
                    unlock(row_locks[row_idx])
                end
            end
        end
    end

    return Z
end

end
