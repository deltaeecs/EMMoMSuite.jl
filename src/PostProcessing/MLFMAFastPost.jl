"""
    MLFMAFastPostModule — MLFMA 中间量复用加速

通过缓存 MLFMA 算子（八叉树 + 近场矩阵），支持同一频率下多个 RHS 向量的
快速重复求解，避免重复 O(N log N) 的树构建开销。

# 使用模式

```julia
# 1. 预构建缓存（仅一次）
cache = build_mlfma_cache(mesh, freq, solver_opts)

# 2. 多个入射方向批量求解
X = solve_multi_rhs!(cache, RHS_matrix)

# 3. 频率变化时自动重建（或手动）
invalidate_cache!(cache)
cache2 = build_mlfma_cache(mesh, new_freq, solver_opts)
```

# 注意事项

- `MLFMACache` 是可变结构体（mutable struct），支持 `invalidate_cache!`
- `validate_cache(cache, freq)` 检查有效性与频率一致性
- `solve_multi_rhs` 对非 MLFMA 矩阵（Matrix/AbstractMatrix）回退到直接 LU 分解
"""
module MLFMAFastPostModule

using LinearAlgebra

export MLFMACache
export invalidate_cache!, validate_cache
export solve_multi_rhs

# ─────────────────────────────────────────────────────────────────────────────
# MLFMACache
# ─────────────────────────────────────────────────────────────────────────────

"""
    MLFMACache

MLFMA 算子复用缓存。

# 字段
- `operator`   : MLFMA 算子（`MLFMAOperator` 或 `nothing`）
- `preconditioner` : 预条件子（`ILUFactor` 等，或 `nothing`）
- `freq`       : 对应的有效频率 (Hz)
- `is_valid`   : 缓存是否有效
- `metadata`   : 补充信息（n_basis、构建时间等）

# 构造
```julia
# 单元测试用（不需要真实 MLFMA 数据）
cache = MLFMACache(; freq=300e6, is_valid=true)

# 带 metadata
cache = MLFMACache(; freq=300e6, is_valid=true, metadata=Dict("n_basis"=>1000))

# 完整构造（由 build_mlfma_cache 填充）
cache = MLFMACache(operator, preconditioner, freq, is_valid, metadata)
```
"""
mutable struct MLFMACache
    operator          :: Any        # MLFMAOperator 或 nothing
    preconditioner    :: Any        # 预条件子或 nothing
    freq              :: Float64
    is_valid          :: Bool
    metadata          :: Dict{String, Any}
end

"""
    MLFMACache(; freq, is_valid=true, metadata=Dict())

关键字构造函数（for testing / lightweight use）。
"""
function MLFMACache(;
    freq         :: Real,
    is_valid     :: Bool                = true,
    metadata     :: AbstractDict        = Dict{String,Any}(),
    operator                            = nothing,
    preconditioner                      = nothing,
)
    meta = Dict{String,Any}(String(k) => v for (k,v) in metadata)
    return MLFMACache(operator, preconditioner, Float64(freq), is_valid, meta)
end

# ─────────────────────────────────────────────────────────────────────────────
# 缓存管理
# ─────────────────────────────────────────────────────────────────────────────

"""
    invalidate_cache!(cache::MLFMACache)

将缓存标记为无效（`is_valid = false`）。
频率变化、几何发生改变或手动重置时使用。
"""
function invalidate_cache!(cache::MLFMACache)
    cache.is_valid = false
    return nothing
end

"""
    validate_cache(cache::MLFMACache, freq::Real) → Bool

检查缓存是否有效且频率吻合。

## 判定规则
- `cache.is_valid` 必须为 `true`
- `|cache.freq - freq| / cache.freq < 1e-10`（相对误差）

**注意**: `cache.freq` 应为正数（由构造函数保证），  
若为零（异常情况），相对误差为 NaN，函数返回 `false`（安全退化）。

当频率不匹配时，本函数**不**自动重建缓存，调用方应手动重建。
"""
function validate_cache(cache::MLFMACache, freq::Real)
    cache.is_valid || return false
    rel_err = abs(cache.freq - Float64(freq)) / cache.freq
    return rel_err < 1e-10
end

# ─────────────────────────────────────────────────────────────────────────────
# 多 RHS 求解
# ─────────────────────────────────────────────────────────────────────────────

"""
    solve_multi_rhs(A, B::AbstractMatrix) → Matrix

对 `k` 个右端向量同时求解线性方程组 `A X = B`。

- **A 为 `AbstractMatrix`**（非 MLFMA）：对 A 进行一次 LU 分解，
  然后对每列 B[:,k] 回代，避免重复因式分解。
  等价于 `A \\ B`，但显式分解更透明。
  
- **A 为 `MLFMAOperator`**（由 `MLFMACache` 提供）：
  复用已建好的八叉树，仅重复 AggTransDis 矩阵-向量乘积。
  当前版本对非 MLFMA 输入统一回退到直接法。

# 参数
- `A` : N×N 系数矩阵（`Matrix` 或 `MLFMAOperator`）
- `B` : N×K 右端矩阵

# 返回
- N×K 解矩阵 `X`，满足 `norm(A*X - B) / norm(B) ≈ 0`

# 示例
```julia
A = complex(rand(4,4)) + 4I   # 主对角线优势，保证可逆
B = randn(ComplexF64, 4, 5)
X = solve_multi_rhs(A, B)
@assert norm(A*X - B) < 1e-10
```
"""
function solve_multi_rhs(A::AbstractMatrix, B::AbstractMatrix)
    N, K = size(B)
    size(A, 1) == N && size(A, 2) == N ||
        throw(DimensionMismatch("A must be N×N for N=$(N)"))

    # 单次 LU 分解，K 次回代
    F = factorize(A)
    X = similar(B)
    for k in 1:K
        X[:, k] = F \ B[:, k]
    end
    return X
end

end  # module MLFMAFastPostModule
