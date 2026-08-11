"""
    FastAlgorithms

快速算法子模块，包含：
- `MLFMA`：多层快速多极算法（含 PMCHW-MLFMA 算子）；
- `Lebedev`：Lebedev 球面积分点集与球面插值（`SHInterp` / `LVI` 等）。

对外导出 `MLFMAOperator`、`get_leaf_intervals`、`PMCHWMLFMAErrorBudget`、
`PMCHWMLFMAOperator`、`assemble_near_field_pmchw` 等。
"""
module FastAlgorithms

include("MLFMA/MLFMA.jl")
include("Lebedev/Lebedev.jl")
include("ACA/ACA.jl")

using .MLFMA
using .Lebedev
using .ACA

export MLFMA, Lebedev, MLFMAOperator, get_leaf_intervals, PMCHWMLFMAErrorBudget, PMCHWMLFMAOperator, assemble_near_field_pmchw
export ACA, LowRankBlock, aca, recompress!, compression_stats

end
