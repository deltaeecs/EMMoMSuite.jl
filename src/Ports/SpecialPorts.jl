"""
    SpecialPorts — DiscretePort 与 CoaxPort

- **DiscretePort**：用于电路节点（SPICE-like 连接），仅定义结构，暂无 MoM 装配。
- **CoaxPort**：同轴探头馈电，TEM 近似特性阻抗 Z_c = (η₀/(2π)) × ln(b/a) / √ε_r。

# 参考
Pozar《Microwave Engineering》 §2.4 同轴传输线。
"""

# ─────────────────────────────────────────────────────────────────────────────
# DiscretePort
# ─────────────────────────────────────────────────────────────────────────────

"""
    DiscretePort{CT<:Complex} <: AbstractPort

离散连接端口（用于电路节点、SPICE 预留接口）。

# 字段
- `id`:            端口编号
- `node_positive`: 正极节点索引（网格节点编号，1-based）
- `node_negative`: 负极节点索引
- `element_type`:  元件类型 `:R | :L | :C | :Z_general`
- `value`:         元件值（Ω | H | F | Ω，取决于 element_type）
"""
struct DiscretePort{CT<:Complex} <: AbstractPort
    id::Int
    node_positive::Int
    node_negative::Int
    element_type::Symbol
    value::CT
end

"""
    DiscretePort(id, node_pos, node_neg, type, value) → DiscretePort

构造离散端口。`value` 支持 Real 或 Complex（自动转为 ComplexF64）。
"""
function DiscretePort(
    id::Int,
    node_positive::Int,
    node_negative::Int,
    element_type::Symbol,
    value::Number,
)
    return DiscretePort{ComplexF64}(id, node_positive, node_negative, element_type, ComplexF64(value))
end

# DiscretePort 暂无电流/电压接口（需要节点电压求解器，超出 Phase 20 范围）
port_current(::DiscretePort, ::AbstractVector) =
    error("DiscretePort port_current requires nodal analysis solver (not yet implemented)")
port_voltage(::DiscretePort, ::AbstractVector) =
    error("DiscretePort port_voltage requires nodal analysis solver (not yet implemented)")

# ─────────────────────────────────────────────────────────────────────────────
# CoaxPort
# ─────────────────────────────────────────────────────────────────────────────

"""
    CoaxPort{FT<:AbstractFloat} <: AbstractPort

同轴探头端口（TEM 近似）。

# 字段
- `id`:             端口编号
- `inner_radius`:   内导体半径 a（m）
- `outer_radius`:   外导体半径 b（m）
- `fill_eps_r`:     填充介质相对介电常数（无损，实数；默认 1.0）

# 特性阻抗公式（TEM 模式）
    Z_c = (η₀ / (2π)) × ln(b/a) / √ε_r

其中 η₀ = 376.730313 Ω 为真空波阻抗。
"""
struct CoaxPort{FT<:AbstractFloat} <: AbstractPort
    id::Int
    inner_radius::FT
    outer_radius::FT
    fill_eps_r::FT
end

"""
    CoaxPort(id, inner_radius, outer_radius, fill_eps_r=1.0) → CoaxPort

构造同轴端口。

# 参数
- `id`:           端口编号
- `inner_radius`: 内导体半径 a（m）
- `outer_radius`: 外导体半径 b（m）
- `fill_eps_r`:   填充介质 ε_r（默认 1.0）
"""
function CoaxPort(
    id::Int,
    inner_radius::Number,
    outer_radius::Number,
    fill_eps_r::Number = 1.0,
)
    FT = promote_type(typeof(float(inner_radius)), typeof(float(outer_radius)), typeof(float(fill_eps_r)))
    return CoaxPort{FT}(id, FT(inner_radius), FT(outer_radius), FT(fill_eps_r))
end

"""
    coax_impedance(port::CoaxPort) → Float64

计算同轴端口特性阻抗（TEM 近似）：

    Z_c = (η₀ / (2π)) × ln(b/a) / √ε_r

# 精度
对于 b/a > 1（正常同轴），误差仅来自浮点舍入（< 1e-10 相对误差）。
"""
function coax_impedance(p::CoaxPort)
    η₀ = 376.730313461  # 真空波阻抗（Ω），与 Constants.eta0 一致
    # Validate coaxial geometry: 0 < inner < outer
    @assert p.inner_radius > 0 "CoaxPort: inner_radius must be > 0"
    @assert p.outer_radius > p.inner_radius "CoaxPort: outer_radius must be > inner_radius"
    @assert p.fill_eps_r > 0 "CoaxPort: fill_eps_r must be > 0"
    return (η₀ / (2π)) * log(p.outer_radius / p.inner_radius) / sqrt(p.fill_eps_r)
end

# CoaxPort 同样不直接对应 MoM 基函数，端口电流/电压预留
port_current(::CoaxPort, ::AbstractVector) =
    error("CoaxPort port_current: need mesh edge association (not yet implemented)")
port_voltage(::CoaxPort, ::AbstractVector) =
    error("CoaxPort port_voltage: need mesh edge association (not yet implemented)")
