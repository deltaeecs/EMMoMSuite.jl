"""
    NearFieldAdvancedModule — 近场结构化输出容器

提供 `NearFieldGrid`（切面网格）和 `NearFieldLine`（切线采样）两种结构体,
可与 VTK 导出函数配合。
"""
module NearFieldAdvancedModule

using StaticArrays
using LinearAlgebra

export NearFieldGrid, NearFieldLine

# ─────────────────────────────────────────────────────────────────────────────
# NearFieldGrid
# ─────────────────────────────────────────────────────────────────────────────

"""
    NearFieldGrid

二维平面切面的近场采样结果。

# 字段
- `coords_u`  : 第一维坐标向量（长度 Nu，单位 m）
- `coords_v`  : 第二维坐标向量（长度 Nv，单位 m）
- `E_field`   : 电场，形状 `(Nu, Nv)`，每元素 `SVector{3, ComplexF64}`
- `H_field`   : 磁场，形状 `(Nu, Nv)`，每元素 `SVector{3, ComplexF64}`
- `freq`      : 频率 (Hz)
- `metadata`  : 补充信息（法向量、偏移量等）

# 构造约束
- `length(coords_u)` 必须等于 `size(E_field, 1)`
- `length(coords_v)` 必须等于 `size(E_field, 2)`
"""
struct NearFieldGrid
    coords_u :: Vector{Float64}
    coords_v :: Vector{Float64}
    E_field  :: Array{SVector{3,ComplexF64}, 2}
    H_field  :: Array{SVector{3,ComplexF64}, 2}
    freq     :: Float64
    metadata :: Dict{String, Any}

    function NearFieldGrid(cu, cv, E, H, freq, meta)
        Nu = length(cu)
        Nv = length(cv)
        size(E) == (Nu, Nv) ||
            throw(ArgumentError(
                "coords_u length=$Nu does not match E_field size $(size(E))"))
        size(H) == (Nu, Nv) ||
            throw(ArgumentError(
                "H_field size $(size(H)) does not match (Nu=$Nu, Nv=$Nv)"))
        return new(
            Float64.(cu), Float64.(cv),   # 强制转换为 Float64
            E, H, Float64(freq), meta
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# NearFieldLine
# ─────────────────────────────────────────────────────────────────────────────

"""
    NearFieldLine

沿任意参数化路径的近场采样结果。

# 字段
- `arc_length` : 弧长参数（单调递增，单位 m）
- `points`     : 三维坐标点
- `E_field`    : 电场矢量数组
- `H_field`    : 磁场矢量数组
- `freq`       : 频率 (Hz)

# 构造约束
- `arc_length`, `points`, `E_field`, `H_field` 长度必须相同
"""
struct NearFieldLine
    arc_length :: Vector{Float64}
    points     :: Vector{SVector{3,Float64}}
    E_field    :: Vector{SVector{3,ComplexF64}}
    H_field    :: Vector{SVector{3,ComplexF64}}
    freq       :: Float64

    function NearFieldLine(arc, pts, E, H, freq)
        N = length(arc)
        length(pts) == N || throw(ArgumentError(
            "points length $(length(pts)) ≠ arc_length length $N"))
        length(E) == N || throw(ArgumentError(
            "E_field length $(length(E)) ≠ arc_length length $N"))
        length(H) == N || throw(ArgumentError(
            "H_field length $(length(H)) ≠ arc_length length $N"))
        return new(Float64.(arc), pts, E, H, Float64(freq))  # arc 强制转换
    end
end

end  # module NearFieldAdvancedModule
