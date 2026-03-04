"""
    RCSBatchModule — RCS 批处理管道

提供 `RCSResult` 容器与便捷查询接口 `rcs_frequency_response` / `rcs_angular_pattern`。
完整批量扫描 API `rcs_batch_sweep` 需要集成 MoM 求解器（见文档）。
"""
module RCSBatchModule

export RCSResult, rcs_frequency_response, rcs_angular_pattern

# ─────────────────────────────────────────────────────────────────────────────
# RCSResult
# ─────────────────────────────────────────────────────────────────────────────

"""
    RCSResult

RCS 批量扫描结果容器。

# 字段
- `freqs`     : 频率向量 (Hz)，长度 Nf
- `theta_inc` : 入射方向极角（rad），长度 Nθ
- `phi_inc`   : 入射方向方位角（rad），长度 Nφ
- `rcs_vv`    : VV 极化 RCS（m²，线性），形状 `(Nf, Nθ, Nφ)`
- `rcs_hh`    : HH 极化 RCS
- `rcs_vh`    : VH（交叉）极化 RCS
- `rcs_hv`    : HV（交叉）极化 RCS

# 构造约束
所有 RCS 数组的形状必须 `== (length(freqs), length(theta_inc), length(phi_inc))`。
"""
struct RCSResult
    freqs     :: Vector{Float64}
    theta_inc :: Vector{Float64}
    phi_inc   :: Vector{Float64}
    rcs_vv    :: Array{Float64, 3}
    rcs_hh    :: Array{Float64, 3}
    rcs_vh    :: Array{Float64, 3}
    rcs_hv    :: Array{Float64, 3}

    function RCSResult(freqs, theta_inc, phi_inc, rcs_vv, rcs_hh, rcs_vh, rcs_hv)
        Nf = length(freqs)
        Nθ = length(theta_inc)
        Nφ = length(phi_inc)
        expected_size = (Nf, Nθ, Nφ)
        for (name, arr) in (
                ("rcs_vv", rcs_vv), ("rcs_hh", rcs_hh),
                ("rcs_vh", rcs_vh), ("rcs_hv", rcs_hv),
        )
            size(arr) == expected_size ||
                throw(ArgumentError(
                    "$name shape $(size(arr)) ≠ expected $(expected_size)"))
        end
        return new(
            collect(Float64, freqs),
            collect(Float64, theta_inc),
            collect(Float64, phi_inc),
            Float64.(rcs_vv), Float64.(rcs_hh),
            Float64.(rcs_vh), Float64.(rcs_hv),
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 内部辅助
# ─────────────────────────────────────────────────────────────────────────────

function _select_rcs_array(result::RCSResult, pol::Symbol)
    pol === :VV && return result.rcs_vv
    pol === :HH && return result.rcs_hh
    pol === :VH && return result.rcs_vh
    pol === :HV && return result.rcs_hv
    throw(ArgumentError("polarization must be :VV, :HH, :VH, or :HV; got $pol"))
end

function _rcs_to_dbsm(rcs_lin::Float64)
    rcs_lin <= 0.0 && return -Inf
    return 10 * log10(rcs_lin)
end

# ─────────────────────────────────────────────────────────────────────────────
# rcs_frequency_response
# ─────────────────────────────────────────────────────────────────────────────

"""
    rcs_frequency_response(result; theta_idx=1, phi_idx=1, polarization=:VV)
        → (freqs, rcs_dBsm)

从 `RCSResult` 提取频率响应曲线。

# 返回
- `freqs`    : 频率向量（Hz）
- `rcs_dBsm` : 对应频点的 RCS（dBsm），长度 Nf

# 参数
- `theta_idx`    : 入射极角索引
- `phi_idx`      : 入射方位角索引
- `polarization` : `:VV`、`:HH`、`:VH`、`:HV`
"""
function rcs_frequency_response(
    result  :: RCSResult;
    theta_idx     :: Int    = 1,
    phi_idx       :: Int    = 1,
    polarization  :: Symbol = :VV,
)
    arr = _select_rcs_array(result, polarization)
    Nf  = size(arr, 1)

    1 <= theta_idx <= size(arr, 2) ||
        throw(ArgumentError("theta_idx=$theta_idx out of range"))
    1 <= phi_idx   <= size(arr, 3) ||
        throw(ArgumentError("phi_idx=$phi_idx out of range"))

    rcs_lin = arr[:, theta_idx, phi_idx]       # Nf
    rcs_dBsm = _rcs_to_dbsm.(rcs_lin)
    return result.freqs, rcs_dBsm
end

# ─────────────────────────────────────────────────────────────────────────────
# rcs_angular_pattern
# ─────────────────────────────────────────────────────────────────────────────

"""
    rcs_angular_pattern(result; freq_idx=1, polarization=:VV, cut_plane=:E)
        → (angles, rcs_dBsm)

从 `RCSResult` 提取角度分布切面。

# 参数
- `freq_idx`    : 频率索引
- `polarization`: `:VV`、`:HH`、`:VH`、`:HV`
- `cut_plane`   : `:E`（沿 θ，phi_idx=1）或 `:H`（沿 φ，theta_idx=1）

# 返回
- `angles`   : 角度向量（rad）
- `rcs_dBsm` : 对应 RCS（dBsm）
"""
function rcs_angular_pattern(
    result       :: RCSResult;
    freq_idx     :: Int    = 1,
    polarization :: Symbol = :VV,
    cut_plane    :: Symbol = :E,
)
    arr = _select_rcs_array(result, polarization)
    Nf  = size(arr, 1)

    1 <= freq_idx <= Nf ||
        throw(ArgumentError("freq_idx=$freq_idx out of range 1:$Nf"))

    if cut_plane === :E
        # E-plane: 沿 θ 切面，固定第一个 φ 角（phi_inc[1]）
        # 注意：这是简化实现，默认主辐射方向在 phi_inc[1] 处，
        # 如需特定 phi 角切面，请直接索引 result.rcs_** 数组
        rcs_lin = arr[freq_idx, :, 1]     # Nθ
        angles  = result.theta_inc
    elseif cut_plane === :H
        # H-plane: 沿 φ 切面，固定第一个 θ 角（theta_inc[1]）
        # 注意：这是简化实现，通常 H-plane 对应最大方向，
        # 如主波束不在 theta_inc[1] 处，请手动选择正确的 theta_idx
        rcs_lin = arr[freq_idx, 1, :]     # Nφ
        angles  = result.phi_inc
    else
        throw(ArgumentError("cut_plane must be :E or :H; got $cut_plane"))
    end

    rcs_dBsm = _rcs_to_dbsm.(rcs_lin)
    return collect(Float64, angles), rcs_dBsm
end

end  # module RCSBatchModule
