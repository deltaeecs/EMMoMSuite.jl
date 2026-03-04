"""
    WavePort — 波端口（矩形波导解析 TE/TM 模式）

当前实现：矩形波导横截面的解析模式（TE_mn / TM_mn），支持：
- 传播常数 β（实数 → 传播；纯虚 → 截止）
- 去嵌入 S 矩阵（相位补偿）

## 物理背景

对于宽 `a`、高 `b` 的矩形波导，TE_mn 模式传播常数：

    β_mn = √(k² − (mπ/a)² − (nπ/b)²)

其中 k = 2πf/c₀。主模 TE10（m=1, n=0）：

    β_10 = √(k² − (π/a)²)

截止频率 f_c = c₀/(2a)（对 TE10）。截止以下 β 为纯虚（指数衰减）。

## 参考
Pozar《Microwave Engineering》§3.3。
"""

# ─────────────────────────────────────────────────────────────────────────────
# PortModes 数据结构
# ─────────────────────────────────────────────────────────────────────────────

"""
    PortModes

矩形波导端口模式求解结果。

# 字段
- `beta`:       传播常数向量（复数，rad/m）；实数 → 传播，纯虚 → 截止
- `kc`:         截止波数向量（rad/m）
- `mode_types`: 模式类型标识（`:TE` 或 `:TM`）
- `m_indices`:  横向指数 m
- `n_indices`:  横向指数 n
- `Z_c`:        模式特性阻抗（Ω）
"""
struct PortModes
    beta::Vector{ComplexF64}
    kc::Vector{Float64}
    mode_types::Vector{Symbol}
    m_indices::Vector{Int}
    n_indices::Vector{Int}
    Z_c::Vector{ComplexF64}
end

# ─────────────────────────────────────────────────────────────────────────────
# WavePort 结构体
# ─────────────────────────────────────────────────────────────────────────────

"""
    WavePort{FT<:AbstractFloat} <: AbstractPort

矩形波导波端口（解析模式）。

# 字段
- `id`:                  端口编号
- `mode`:                激励模式（`:TE10 | :TE20 | :TE01 | :TM11 | ...`）
- `a`:                   波导宽度（m）
- `b`:                   波导高度（m）
- `de_embed_length`:     去嵌入距离 L（m）；0.0 = 无去嵌入
- `reference_impedance`:参考阻抗（Ω，用于 S 参数）
"""
struct WavePort{FT<:AbstractFloat} <: AbstractPort
    id::Int
    mode::Symbol
    a::FT
    b::FT
    de_embed_length::FT
    reference_impedance::FT
end

"""
    WavePort(id; mode, a, b, de_embed_length=0.0, reference_impedance=50.0) → WavePort

构造矩形波导端口。

# 关键字参数
- `mode`:  模式符号，如 `:TE10`（默认）
- `a`:     波导宽度（m）
- `b`:     波导高度（m，默认 a/2）
- `de_embed_length`: 去嵌入距离（m，默认 0.0）
- `reference_impedance`: 参考阻抗（Ω，默认 50.0）
"""
function WavePort(
    id::Int;
    mode::Symbol        = :TE10,
    a::Real             = 1.0,
    b::Real             = a / 2,
    de_embed_length::Real       = 0.0,
    reference_impedance::Real   = 50.0,
)
    FT = promote_type(typeof(float(a)), typeof(float(b)),
                      typeof(float(de_embed_length)), typeof(float(reference_impedance)))
    return WavePort{FT}(id, mode, FT(a), FT(b), FT(de_embed_length), FT(reference_impedance))
end

# WavePort 本身不直接对应单个 RWG 基函数，留待 meshed-port 扩展
port_current(::WavePort, ::AbstractVector) =
    error("WavePort port_current: requires meshed port DOF association (Phase 20 扩展)")
port_voltage(::WavePort, ::AbstractVector) =
    error("WavePort port_voltage: requires meshed port DOF association (Phase 20 扩展)")

# ─────────────────────────────────────────────────────────────────────────────
# 模式参数解析
# ─────────────────────────────────────────────────────────────────────────────

"""
    _parse_mode_indices(mode::Symbol) → (Symbol, Int, Int)

从模式符号提取 (type, m, n)。

# 支持格式
- `:TE10` → (:TE, 1, 0)
- `:TM21` → (:TM, 2, 1)

**限制**：仅支持单位数指数（m, n ∈ 0–9）。
对于 TE110 等双位数指数模式，需扩展解析逻辑。
"""
function _parse_mode_indices(mode::Symbol)
    s = string(mode)
    @assert length(s) >= 4 "模式符号格式应为 :TEmn 或 :TMmn，例如 :TE10"
    type_str = s[1:2]
    m = parse(Int, string(s[3]))
    n = parse(Int, string(s[4]))
    type_sym = type_str == "TE" ? :TE : :TM
    return type_sym, m, n
end

# ─────────────────────────────────────────────────────────────────────────────
# 模式求解（解析）
# ─────────────────────────────────────────────────────────────────────────────

"""
    compute_port_modes(port::WavePort, freq::Real, n_modes::Int=1) → PortModes

计算矩形波导端口的解析模式参数。

对于主模（n_modes=1），仅计算 `port.mode` 对应的单个模式。
对于 n_modes>1，按截止波数升序返回前 n_modes 个 TE 模式（忽略 TM，为简化实现）。

# 精度
使用 Float64 双精度算术，解析公式无数值误差（仅浮点舍入 ~1e-15 量级）。

# 返回
`PortModes` 包含 `n_modes` 个模式的完整参数。
"""
function compute_port_modes(port::WavePort, freq::Real, n_modes::Int = 1)
    c0 = 299792458.0      # m/s
    η₀ = 376.730313461   # 真空波阻抗
    k  = ComplexF64(2π * freq / c0)

    if n_modes == 1
        # 仅计算指定模式
        type_sym, m, n = _parse_mode_indices(port.mode)
        kc   = sqrt((m * π / port.a)^2 + (n * π / port.b)^2)
        beta = sqrt(k^2 - kc^2)
        # 确保 beta 具有非负实部（传播方向约定）
        if real(beta) < 0
            beta = -beta
        end

        # 特性阻抗
        Z_c = if type_sym === :TE
            k * η₀ / beta     # Z_TE = η₀ k / β
        else
            beta * η₀ / k     # Z_TM = η₀ β / k
        end

        return PortModes(
            [beta],
            [kc],
            [type_sym],
            [m],
            [n],
            [Z_c],
        )
    else
        # 计算前 n_modes 个 TE 模式（m=0,1,2...  n=0,1,2...）
        candidates = Tuple{Float64, Int, Int}[]
        for mm in 0:10, nn in 0:10
            mm == 0 && nn == 0 && continue  # TE00 不存在
            kc_mn = sqrt((mm * π / port.a)^2 + (nn * π / port.b)^2)
            push!(candidates, (kc_mn, mm, nn))
        end
        sort!(candidates, by = x -> x[1])

        betas  = ComplexF64[]
        kcs    = Float64[]
        types  = Symbol[]
        ms     = Int[]
        ns     = Int[]
        Zcs    = ComplexF64[]

        for (kc_mn, mm, nn) in candidates[1:min(n_modes, length(candidates))]
            beta_mn = sqrt(k^2 - kc_mn^2)
            real(beta_mn) < 0 && (beta_mn = -beta_mn)
            Z_c_mn = k * η₀ / beta_mn   # TE 模式
            push!(betas, beta_mn)
            push!(kcs,   kc_mn)
            push!(types, :TE)
            push!(ms,    mm)
            push!(ns,    nn)
            push!(Zcs,   Z_c_mn)
        end

        return PortModes(betas, kcs, types, ms, ns, Zcs)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 去嵌入
# ─────────────────────────────────────────────────────────────────────────────

"""
    de_embed_s_matrix(S_raw, ports::Vector{WavePort}, freq::Real) → Matrix{ComplexF64}

将原始 S 矩阵按各端口的 `de_embed_length` 进行相位去嵌入。

去嵌入公式：
    S_ij_deembed = S_ij_raw × exp(j β_i L_i) × exp(j β_j L_j)

其中 L_i 为端口 i 的去嵌入距离，β_i 为端口 i 的主模传播常数。

# 参数
- `S_raw`:  原始 S 矩阵（含传输线相移）
- `ports`:  端口列表（需与 S_raw 行/列对应）
- `freq`:   频率（Hz）
"""
function de_embed_s_matrix(
    S_raw::AbstractMatrix,
    ports::Vector{<:WavePort},
    freq::Real,
)
    N = length(ports)
    @assert size(S_raw) == (N, N) "S_raw 维度须与端口数一致"

    # 计算每个端口的去嵌入相位因子 exp(j β L)
    phase = ComplexF64[]
    for p in ports
        modes = compute_port_modes(p, freq, 1)
        β  = modes.beta[1]
        push!(phase, exp(im * β * p.de_embed_length))
    end

    # 应用去嵌入：S_deembed[i,j] = S_raw[i,j] * phase[i] * phase[j]
    S_deembed = copy(ComplexF64.(S_raw))
    for j in 1:N, i in 1:N
        S_deembed[i, j] *= phase[i] * phase[j]
    end
    return S_deembed
end
