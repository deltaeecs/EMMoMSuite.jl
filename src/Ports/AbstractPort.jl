"""
    AbstractPort — 端口体系抽象接口

所有端口类型的公共基类型，及其通用接口方法。

# 接口要求（子类型必须定义）
- `port_id(port)` → Int
- `port_current(port, I_coeff)` → ComplexF64
- `port_voltage(port, I_coeff)` → ComplexF64

# 可选实现（有默认基于 V/I 的版本）
- `port_power(port, I_coeff)` → Float64
"""

"""
    AbstractPort

所有端口（集总端口、波端口、同轴端口等）的公共抽象基类型。
"""
abstract type AbstractPort end

"""
    port_id(port) → Int

返回端口编号（唯一标识符）。

**约定**：默认实现直接读取 `port.id` 字段。
所有 `AbstractPort` 子类型**必须**定义名为 `id::Int` 的字段，
或覆盖本方法提供自定义实现。
"""
port_id(p::AbstractPort) = p.id

"""
    port_current(port, I_coeff) → ComplexF64

从阻抗矩阵求解的基函数系数向量 `I_coeff` 中提取端口电流。
"""
function port_current(::AbstractPort, ::AbstractVector)
    error("port_current not implemented for this port type")
end

"""
    port_voltage(port, I_coeff) → ComplexF64

从基函数系数向量 `I_coeff` 提取端口电压（= Z_load × I_port）。
"""
function port_voltage(::AbstractPort, ::AbstractVector)
    error("port_voltage not implemented for this port type")
end

"""
    port_power(port, I_coeff) → Float64

计算端口吸收功率（时均值）：P = ½ Re(V_port · I_port*)。
默认实现基于 `port_voltage` 和 `port_current`。
"""
function port_power(p::AbstractPort, I_coeff::AbstractVector)
    V = port_voltage(p, I_coeff)
    I = port_current(p, I_coeff)
    return 0.5 * real(V * conj(I))
end
