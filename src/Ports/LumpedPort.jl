"""
    LumpedPort — 集总端口

在 RWG MoM 框架中，集总端口对应于网格中某条边（edge）处插入集总阻抗/电源。

# 物理模型
- 在阻抗矩阵的 `(edge_idx, edge_idx)` 对角元素添加集总阻抗 Z_load。
- 对于电压源：在激励向量的 `edge_idx` 位置添加 V_source。
- 端口电流 I_port = I_coeff[edge_idx]（基函数系数）。
- 端口电压 V_port = Z_load × I_port。

# 参考
Harrington《Field Computation by Moment Methods》§4.3 集总负载。
"""

"""
    LumpedPort{CT<:Complex} <: AbstractPort

集总端口。

# 字段
- `id`:               端口唯一编号（正整数）
- `edge_idx`:         对应 RWG 基函数/边索引（1-based）
- `impedance`:        端口阻抗 Z_load（Ω，复数，支持 R/L/C）
- `port_type`:        `:voltage_source | :current_source | :load`
- `voltage_amplitude`:激励电压幅度（电压源类型）
- `current_amplitude`:激励电流幅度（电流源类型）
"""
struct LumpedPort{CT<:Complex} <: AbstractPort
    id::Int
    edge_idx::Int
    impedance::CT
    port_type::Symbol
    voltage_amplitude::CT
    current_amplitude::CT
end

"""
    lumped_port(id, edge_idx; impedance, type, voltage, current) → LumpedPort

便捷构造函数。

# 参数
- `id`:       端口编号
- `edge_idx`: 边索引（1-based）
- `impedance`:端口阻抗（默认 50.0 Ω）
- `type`:     `:voltage_source | :current_source | :load`（默认 `:load`）
- `voltage`:  激励电压（默认 0.0）
- `current`:  激励电流（默认 0.0）
"""
function lumped_port(
    id::Int,
    edge_idx::Int;
    impedance::Number   = 50.0,
    type::Symbol        = :load,
    voltage::Number     = 0.0,
    current::Number     = 0.0,
)
    CT = ComplexF64
    return LumpedPort{CT}(id, edge_idx, CT(impedance), type, CT(voltage), CT(current))
end

# ─────────────────────────────────────────────────────────────────────────────
# 接口实现
# ─────────────────────────────────────────────────────────────────────────────

"""
    port_current(port::LumpedPort, I_coeff) → ComplexF64

端口电流 = MoM 基函数系数 I_coeff[edge_idx]（直接对应，无额外缩放）。
"""
port_current(p::LumpedPort, I_coeff::AbstractVector) = I_coeff[p.edge_idx]

"""
    port_voltage(port::LumpedPort, I_coeff) → ComplexF64

端口电压 = Z_load × I_port。
"""
port_voltage(p::LumpedPort, I_coeff::AbstractVector) = p.impedance * I_coeff[p.edge_idx]

# ─────────────────────────────────────────────────────────────────────────────
# 阻抗矩阵与激励向量贡献
# ─────────────────────────────────────────────────────────────────────────────

"""
    assemble_lumped_port_impedance!(Z, port::LumpedPort)

将集总端口的阻抗贡献加到阻抗矩阵 `Z` 的对角元素。

# 物理意义
对角元素 Z[edge_idx, edge_idx] += Z_load：
将集总阻抗串联到基函数 `edge_idx` 所在的路径上。
"""
function assemble_lumped_port_impedance!(Z::AbstractMatrix, p::LumpedPort)
    Z[p.edge_idx, p.edge_idx] += p.impedance
    return nothing
end

"""
    add_port_excitation!(V, port::LumpedPort)

将集总端口的激励贡献加到右端向量 `V`。

- `:voltage_source`：V[edge_idx] += voltage_amplitude
- `:current_source`：（当前不支持，保留扩展）
- `:load`：不添加激励（纯负载，无源）
"""
function add_port_excitation!(V::AbstractVector, p::LumpedPort)
    if p.port_type === :voltage_source
        V[p.edge_idx] += p.voltage_amplitude
    elseif p.port_type === :current_source
        # 电流源：J = I/Z（诺顿等效变换），预留接口
        V[p.edge_idx] += p.current_amplitude * p.impedance
    end
    # :load 类型无激励贡献
    return nothing
end
