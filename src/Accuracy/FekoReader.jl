"""
    FekoReader — 解析 Feko 远场 RCS CSV 文件

# 文件格式说明（MoM_AllinOne deps/compare_feko/*.csv）

每个文件包含 1443 行：
- 第 1 行: 标题（以 "THETA" 开头，115 字符）
- 第 2–1443 行: 数据行（每行 ~113 字符，包含所有字段）

每行数据字段（空格分隔）：
  1: THETA (deg)
  2: PHI (deg)
  3: E_theta 幅值
  4: E_theta 相位 (deg)
  5: E_phi 幅值
  6: E_phi 相位 (deg)
  7: RCS (m²)    ← 我们使用此列
  8: 轴比
  9: 倾斜角
 10: 极化类型

每个文件含两个 φ 切面（φ=0° 和 φ=90°），每切面 721 个 θ 点（-180° 到 180°，步长 0.5°）。
"""
module FekoReader

export read_feko_rcs, split_phi_cuts

"""
    read_feko_rcs(filepath) -> (theta_deg, phi_deg, rcs_sqm, rcs_dBsm)

解析 Feko 远场 RCS CSV 文件。

# 返回
- `theta_deg  :: Vector{Float64}` — θ 角度（度），长度 1442
- `phi_deg    :: Vector{Float64}` — φ 角度（度），长度 1442
- `rcs_sqm    :: Vector{Float64}` — RCS（m²），长度 1442
- `rcs_dBsm   :: Vector{Float64}` — RCS（dBsm = 10·log₁₀(m²)），长度 1442

# 异常
- `ArgumentError`: 文件不存在
- `ErrorException`: 文件格式不符（无法解析）
"""
function read_feko_rcs(filepath::AbstractString)
    isfile(filepath) || throw(ArgumentError("Feko 文件不存在: $filepath"))

    theta_deg  = Float64[]
    phi_deg    = Float64[]
    rcs_sqm    = Float64[]

    open(filepath, "r") do io
        for line in eachline(io)
            stripped = strip(line)
            isempty(stripped) && continue

            # 跳过标题行：包含 "THETA" 或 "axial"
            (startswith(stripped, "THETA") || startswith(stripped, "axial")) && continue

            parts = split(stripped)
            length(parts) < 7 && continue

            # 尝试解析为浮点，跳过非数字行（极化标签如 "LINEAR"）
            θ = tryparse(Float64, parts[1])
            isnothing(θ) && continue
            φ = tryparse(Float64, parts[2])
            isnothing(φ) && continue
            r = tryparse(Float64, parts[7])
            isnothing(r) && continue

            push!(theta_deg, θ)
            push!(phi_deg,   φ)
            push!(rcs_sqm,   r)
        end
    end

    isempty(rcs_sqm) && error("未能从文件解析任何数据: $filepath")

    rcs_dBsm = 10.0 .* log10.(rcs_sqm)

    return theta_deg, phi_deg, rcs_sqm, rcs_dBsm
end

"""
    split_phi_cuts(theta_deg, phi_deg, rcs_dBsm; tol=0.1)
        -> Dict{Float64, NamedTuple}

按 φ 值分组，返回 `φ → (theta, rcs_dBsm)` 的有序字典。
每个切面的 θ 值按升序排列。

# 参数
- `tol`: φ 值分组容差（度），默认 0.1°

# 返回
`Dict{Float64, NamedTuple{(:theta, :rcs_dBsm), Tuple{Vector{Float64}, Vector{Float64}}}}`
键为四舍五入到 0.01° 的 φ 值。
"""
function split_phi_cuts(
    theta_deg::Vector{Float64},
    phi_deg::Vector{Float64},
    rcs_dBsm::Vector{Float64};
    tol::Float64 = 0.1,
)
    # 找到所有不同的 φ 切面值（四舍五入到 0.1 度）
    phi_rounded = round.(phi_deg, digits = 1)
    phi_unique  = sort(unique(phi_rounded))

    result = Dict{Float64,NamedTuple{(:theta, :rcs_dBsm),Tuple{Vector{Float64},Vector{Float64}}}}()

    for φ_val in phi_unique
        mask = abs.(phi_rounded .- φ_val) .< tol
        th   = theta_deg[mask]
        rcs  = rcs_dBsm[mask]

        # 按 θ 升序排列
        order = sortperm(th)
        result[φ_val] = (theta = th[order], rcs_dBsm = rcs[order])
    end

    return result
end

end # module FekoReader
