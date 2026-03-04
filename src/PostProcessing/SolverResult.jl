"""
    SolverResultModule — MoM 求解结果封装

将 `(basis, I_coeffs, freq)` 三元组包装为具名结构体，同时携带可选 metadata。

提供向后兼容辅助，使现有 `(basis, I_coeffs)` 风格调用过渡顺畅。
"""
module SolverResultModule

export SolverResult

"""
    SolverResult{BT, CT <: Number}

MoM 求解结果容器。

# 字段
- `basis`      : 基函数对象（RWGBasis、SWGBasis 等；单元测试可传 `nothing`）
- `I_coeffs`   : 电流系数向量（复数）
- `freq`       : 频率 (Hz)，必须 > 0
- `metadata`   : 可选补充信息（方程类型、网格尺寸等）

# 构造
```julia
sr = SolverResult(basis, I, 1e9)
sr = SolverResult(basis, I, 1e9; metadata=Dict("equation"=>"EFIE"))
```
"""
struct SolverResult{BT, CT <: Number}
    basis    :: BT
    I_coeffs :: Vector{CT}
    freq     :: Float64
    metadata :: Dict{String, Any}
end

"""
    SolverResult(basis, I_coeffs, freq; metadata=Dict{String,Any}())

外部构造函数——验证 freq > 0。
"""
function SolverResult(
    basis,
    I_coeffs::Vector{CT},
    freq::Real;
    metadata::AbstractDict = Dict{String,Any}(),
) where {CT <: Number}
    freq > 0 || throw(ArgumentError("freq must be positive, got $freq"))
    meta = Dict{String,Any}(String(k) => v for (k,v) in metadata)
    return SolverResult{typeof(basis), CT}(basis, I_coeffs, Float64(freq), meta)
end

# ── 便捷接口 ──────────────────────────────────────────────────────────────────

"""
    Base.length(sr::SolverResult) → Int

返回系数数量（等同 `num_basis(sr.basis)`）。
"""
Base.length(sr::SolverResult) = length(sr.I_coeffs)

end  # module SolverResultModule
