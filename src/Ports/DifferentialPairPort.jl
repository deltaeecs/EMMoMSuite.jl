"""
    DifferentialPairPort — 差分端口与混合模 S 参数转换

将单端（single-ended）N 端口 S 矩阵转换为混合模（mixed-mode）S 矩阵。

# 混合模变换
对于 N/2 个差分对（N 必须为偶数），变换矩阵：

    M = (1/√2) × blockdiag([1 1; 1 -1], ...)

其中每个 2×2 块对应一个差分对（正端口 → 行 2i-1，负端口 → 行 2i）。

转换公式：
    S_mm = M × S_se × M†

分量：
- S_dd（差分-差分）：远场信号质量
- S_cc（共模-共模）：EMI 指标
- S_cd（转换分量）：CMRR 来源（理想情况应为 0）

# 参考
Pozar《Microwave Engineering》§4.8；Bogatin《Signal Integrity》§9.
"""

# ─────────────────────────────────────────────────────────────────────────────
# 结构体
# ─────────────────────────────────────────────────────────────────────────────

"""
    DifferentialPairPort <: AbstractPort

差分端口（关联正负两个 LumpedPort 或逻辑端口）。

# 字段
- `id`:            差分对编号
- `port_positive`: 正线端口 id（关联已有 LumpedPort/DiscretePort）
- `port_negative`: 负线端口 id
"""
struct DifferentialPairPort <: AbstractPort
    id::Int
    port_positive::Int
    port_negative::Int
end

# DifferentialPairPort 本身无直接 MoM 电流/电压接口
port_current(::DifferentialPairPort, ::AbstractVector) =
    error("DifferentialPairPort: use underlying positive/negative ports for current extraction")
port_voltage(::DifferentialPairPort, ::AbstractVector) =
    error("DifferentialPairPort: use underlying positive/negative ports for voltage extraction")

# ─────────────────────────────────────────────────────────────────────────────
# 混合模变换矩阵
# ─────────────────────────────────────────────────────────────────────────────

"""
    mixed_mode_transform_matrix(pairs::Vector{DifferentialPairPort}) → Matrix{ComplexF64}

构造 N×N 混合模变换矩阵 M（N = 2 × 端口对数）。

# 矩阵结构
对于 K 个差分对（共 2K 端口），M 为 2K×2K 分块对角矩阵，
每个 2×2 块：

    M_k = (1/√2) × [1  1]
                    [1 -1]

满足酉性：M · M† = I。

# 端口排列约束
调用前，必须保证 `pairs` 中 port_positive < port_negative，
且端口编号从 1 到 2K 连续（与 S_se 的列/行顺序一致）。
"""
function mixed_mode_transform_matrix(pairs::Vector{DifferentialPairPort})
    K = length(pairs)
    N = 2K
    M = zeros(ComplexF64, N, N)
    inv_sqrt2 = 1.0 / sqrt(2.0)
    for k in 1:K
        # 差分模（第 2k-1 行）对应正极
        i_d = 2k - 1
        # 共模（第 2k 行）对应正极与负极之和
        i_c = 2k
        # 正极端口 → 列 2k-1；负极端口 → 列 2k
        M[i_d, 2k-1] =  inv_sqrt2
        M[i_d, 2k  ] = -inv_sqrt2   # 差分：V_d = (V+ - V−)/√2
        M[i_c, 2k-1] =  inv_sqrt2
        M[i_c, 2k  ] =  inv_sqrt2   # 共模：V_c = (V+ + V−)/√2
    end
    return M
end

# ─────────────────────────────────────────────────────────────────────────────
# 混合模 S 矩阵转换
# ─────────────────────────────────────────────────────────────────────────────

"""
    convert_to_mixed_mode(S_se, pairs) → Matrix{ComplexF64}

将单端口 S 矩阵 `S_se`（2K×2K）转换为混合模 S 矩阵。

    S_mm = M × S_se × M†

# 参数
- `S_se`:  单端 S 矩阵（2K×2K，ComplexF64）
- `pairs`: 差分对向量（K 个 DifferentialPairPort）
"""
function convert_to_mixed_mode(S_se::AbstractMatrix, pairs::Vector{DifferentialPairPort})
    M = mixed_mode_transform_matrix(pairs)
    return M * S_se * M'
end

# ─────────────────────────────────────────────────────────────────────────────
# 分量提取
# ─────────────────────────────────────────────────────────────────────────────

"""
    sdd_matrix(S_mm, pairs) → Matrix{ComplexF64}

从混合模 S 矩阵提取差分-差分子块（行/列奇数部分）。
"""
function sdd_matrix(S_mm::AbstractMatrix, pairs::Vector{DifferentialPairPort})
    K = length(pairs)
    idx = [2k-1 for k in 1:K]
    return S_mm[idx, idx]
end

"""
    scc_matrix(S_mm, pairs) → Matrix{ComplexF64}

从混合模 S 矩阵提取共模-共模子块（行/列偶数部分）。
"""
function scc_matrix(S_mm::AbstractMatrix, pairs::Vector{DifferentialPairPort})
    K = length(pairs)
    idx = [2k for k in 1:K]
    return S_mm[idx, idx]
end

"""
    scd_matrix(S_mm, pairs) → Matrix{ComplexF64}

从混合模 S 矩阵提取差分-共模转换子块（CMRR 来源）。
理想对称传输线时应全为零。
"""
function scd_matrix(S_mm::AbstractMatrix, pairs::Vector{DifferentialPairPort})
    K = length(pairs)
    d_idx = [2k-1 for k in 1:K]
    c_idx = [2k   for k in 1:K]
    return S_mm[d_idx, c_idx]
end
