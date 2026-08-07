"""
    Ports — EMMoMSuite 端口体系

Phase 20 实现的端口类型：
- `LumpedPort`:            集总端口（R/L/C/电压源/电流源）
- `DiscretePort`:          离散节点端口（SPICE 预留接口）
- `CoaxPort`:              同轴端口（TEM 模式，解析 Z_c）
- `DifferentialPairPort`:  差分对（混合模 S 参数转换）
- `WavePort`:              波端口（矩形波导解析 TE/TM 模式）

通用工具：
- `assemble_lumped_port_impedance!`: 阻抗矩阵贡献
- `add_port_excitation!`:            激励向量贡献
- `port_voltage/current/power`:      端口电量提取
- `extract_s_matrix`:                Z→S 矩阵转换
- `SParameterData`:                  多端口 S 参数数据结构
- `write_touchstone / read_touchstone`: Touchstone 1.0 I/O
- `coax_impedance`:                  同轴特性阻抗
- `mixed_mode_transform_matrix`:     混合模变换矩阵
- `convert_to_mixed_mode`:           单端→混合模转换
- `sdd_matrix / scc_matrix / scd_matrix`: 混合模分量提取
- `compute_port_modes`:              矩形波导模式参数
- `de_embed_s_matrix`:               去嵌入 S 矩阵
"""
module PortsModule

using LinearAlgebra
using ..CoreModule: Constants

include("AbstractPort.jl")
include("LumpedPort.jl")
include("SpecialPorts.jl")
include("DifferentialPairPort.jl")
include("SParameterExtraction.jl")
include("WavePort.jl")

export AbstractPort, port_id, port_current, port_voltage, port_power

# LumpedPort
export LumpedPort, lumped_port
export assemble_lumped_port_impedance!, add_port_excitation!

# DiscretePort + CoaxPort
export DiscretePort, CoaxPort, coax_impedance

# DifferentialPairPort
export DifferentialPairPort
export mixed_mode_transform_matrix, convert_to_mixed_mode
export sdd_matrix, scc_matrix, scd_matrix

# SParameterExtraction
export SParameterData, sparam_data, extract_s_matrix
export write_touchstone, read_touchstone

# WavePort
export WavePort, PortModes
export compute_port_modes, de_embed_s_matrix

end # module PortsModule
