"""
    ReferenceData — Phase 14 精度对比参考数据

提供以下解析参考：
1. Mie PEC 球 RCS（调用 EMSuite 内置 MieSeries）
2. Mie 均匀介质球 RCS（调用 EMSuite 内置 calculate_mie_rcs_dielectric_sphere）
3. 半波偶极子解析输入阻抗（Balanis §8.1）
4. 谐振偶极子（0.47λ）解析阻抗
5. 半波偶极子远场方向图
"""
module ReferenceData

using LinearAlgebra
using Statistics: mean

# 使用 EMSuite Utilities 中已有的 Mie 级数
import ..Accuracy  # 使用父包的 Utilities
import ...Geometry: read_nas_mesh
# 直接调用 EMSuite 模块级导出（在运行时可访问）

export mie_pec_rcs_dBsm, mie_pec_bistatic_rcs_dBsm,
       mie_dielectric_rcs_dBsm, mie_dielectric_bistatic_rcs_dBsm,
       dipole_halfwave_Zin_analytic, dipole_resonant_Zin_analytic,
       dipole_halfwave_farfield_analytic,
       extract_sphere_radius

# ─────────────────────────────────────────────────────────────────────────────
# 辅助：从 Nastran 网格文件提取球半径
# ─────────────────────────────────────────────────────────────────────────────

"""
    extract_sphere_radius(nas_file) -> radius_m

读取 Nastran GRID 节点坐标，计算节点到原点的均值距离作为球半径。
"""
function extract_sphere_radius(nas_file::AbstractString)
    isfile(nas_file) || throw(ArgumentError("网格文件不存在: $nas_file"))

    # Prefer the project's Nastran reader since it already handles GRID*,
    # fixed/free field and compressed scientific notation robustly.
    try
        mesh = read_nas_mesh(String(nas_file); scale = 1.0)
        nodes = getfield(mesh, :node)
        if size(nodes, 1) == 3 && size(nodes, 2) > 0
            rs = vec(sqrt.(sum(abs2, nodes; dims = 1)))
            return mean(rs)
        end
    catch
        # Fallback to lightweight text parsing below.
    end

    rs = Float64[]
    open(nas_file, "r") do io
        for line in eachline(io)
            sline = strip(line)
            startswith(sline, "GRID") || continue
            # Nastran 格式：GRID    id    cp    x     y     z
            parts = split(sline)
            length(parts) < 6 && continue
            x = tryparse(Float64, parts[end-2])
            y = tryparse(Float64, parts[end-1])
            z = tryparse(Float64, parts[end])
            (isnothing(x) || isnothing(y) || isnothing(z)) && continue
            push!(rs, sqrt(x^2 + y^2 + z^2))
        end
    end
    isempty(rs) && error("未能从网格文件提取节点坐标: $nas_file")
    return mean(rs)
end

# ─────────────────────────────────────────────────────────────────────────────
# Mie PEC 球
# ─────────────────────────────────────────────────────────────────────────────

"""
    mie_pec_rcs_dBsm(radius_m, freq_hz, theta_rad_vec) -> rcs_dBsm

计算 PEC 球 Mie 级数 RCS（E 面 S₂ 分量），返回 dBsm。

# 参数
- `theta_rad_vec`: 散射角（弧度），0=前向，π=后向
"""
function mie_pec_rcs_dBsm(radius_m::Real, freq_hz::Real, theta_rad_vec::AbstractVector)
    # 动态调用 EMSuite 顶层导出函数（运行时加载，避免循环依赖）
    # calculate_mie_rcs_pec_sphere 定义在 EMSuite.Utilities.MieSeries
    rcs_sqm = _call_mie_pec(radius_m, freq_hz, theta_rad_vec)
    return 10.0 .* log10.(max.(rcs_sqm, 1e-100))
end

"""
    mie_pec_rcs_dBsm(mesh_file, freq_hz, theta_rad_vec) -> rcs_dBsm

从网格文件提取球半径后计算。
"""
function mie_pec_rcs_dBsm(mesh_file::AbstractString, freq_hz::Real, theta_rad_vec::AbstractVector)
    r = extract_sphere_radius(mesh_file)
    return mie_pec_rcs_dBsm(r, freq_hz, theta_rad_vec)
end

"""
    mie_pec_bistatic_rcs_dBsm(radius_m, freq_hz, theta_obs, phi_obs, theta_inc, phi_inc, polarization)
        -> rcs_dBsm

将 EMSuite 使用的全局观测角 `(theta_obs, phi_obs)` 映射为相对入射方向的散射角，
并使用 PEC 球 Mie 全极化解重建总双站 RCS。

# 参数
- `theta_obs`: 全局观测极角（弧度）
- `phi_obs`: 全局观测方位角（弧度，标量或与 `theta_obs` 等长向量）
- `theta_inc`, `phi_inc`: 入射方向的全局球坐标（弧度）
- `polarization`: 入射电场极化向量（自动投影到垂直于传播方向的平面）
"""
function mie_pec_bistatic_rcs_dBsm(
    radius_m::Real,
    freq_hz::Real,
    theta_obs::AbstractVector,
    phi_obs,
    theta_inc::Real,
    phi_inc::Real,
    polarization::AbstractVector,
)
    phi_vec = _expand_phi_obs(phi_obs, length(theta_obs))
    theta_scat, phi_local = _bistatic_angles(theta_obs, phi_vec, theta_inc, phi_inc, polarization)
    rcs_s2, rcs_s1 = _call_mie_pec_fullpol(radius_m, freq_hz, theta_scat)

    rcs_total = similar(rcs_s2)
    @inbounds for i in eachindex(rcs_total)
        cφ = cos(phi_local[i])
        sφ = sin(phi_local[i])
        rcs_total[i] = rcs_s2[i] * cφ^2 + rcs_s1[i] * sφ^2
    end
    return 10.0 .* log10.(max.(rcs_total, 1e-100))
end

function mie_pec_bistatic_rcs_dBsm(
    mesh_file::AbstractString,
    freq_hz::Real,
    theta_obs::AbstractVector,
    phi_obs,
    theta_inc::Real,
    phi_inc::Real,
    polarization::AbstractVector,
)
    r = extract_sphere_radius(mesh_file)
    return mie_pec_bistatic_rcs_dBsm(r, freq_hz, theta_obs, phi_obs, theta_inc, phi_inc, polarization)
end

# ─────────────────────────────────────────────────────────────────────────────
# Mie 均匀介质球
# ─────────────────────────────────────────────────────────────────────────────

"""
    mie_dielectric_rcs_dBsm(radius_m, freq_hz, eps_r, mu_r, theta_rad_vec) -> rcs_dBsm

计算均匀介质球 Mie 级数的非偏振总 RCS，返回 dBsm。

# 参数
- `eps_r`: 相对介电常数（实数=无损，复数=有损）
- `mu_r`:  相对磁导率（默认 1.0）
- `theta_rad_vec`: 散射角（弧度），0=前向，π=后向
"""
function mie_dielectric_rcs_dBsm(
    radius_m::Real,
    freq_hz::Real,
    eps_r::Number,
    mu_r::Number,
    theta_rad_vec::AbstractVector,
)
    _, _, rcs_unpol = _call_mie_dielectric_fullpol(radius_m, freq_hz, theta_rad_vec, eps_r, mu_r)
    return 10.0 .* log10.(max.(rcs_unpol, 1e-100))
end

"""
    mie_dielectric_bistatic_rcs_dBsm(radius_m, freq_hz, eps_r, mu_r, theta_obs, phi_obs, theta_inc, phi_inc, polarization)
        -> rcs_dBsm

将 EMSuite 使用的全局观测角 `(theta_obs, phi_obs)` 映射为相对入射方向的散射角，
并使用均匀介质球 Mie 全极化解重建总双站 RCS。
"""
function mie_dielectric_bistatic_rcs_dBsm(
    radius_m::Real,
    freq_hz::Real,
    eps_r::Number,
    mu_r::Number,
    theta_obs::AbstractVector,
    phi_obs,
    theta_inc::Real,
    phi_inc::Real,
    polarization::AbstractVector,
)
    phi_vec = _expand_phi_obs(phi_obs, length(theta_obs))
    theta_scat, phi_local = _bistatic_angles(theta_obs, phi_vec, theta_inc, phi_inc, polarization)
    rcs_s2, rcs_s1, _ = _call_mie_dielectric_fullpol(radius_m, freq_hz, theta_scat, eps_r, mu_r)

    rcs_total = similar(rcs_s2)
    @inbounds for i in eachindex(rcs_total)
        cφ = cos(phi_local[i])
        sφ = sin(phi_local[i])
        rcs_total[i] = rcs_s2[i] * cφ^2 + rcs_s1[i] * sφ^2
    end
    return 10.0 .* log10.(max.(rcs_total, 1e-100))
end

# ─────────────────────────────────────────────────────────────────────────────
# 内部辅助：延迟调用 EMSuite 顶层 Mie 函数
# ─────────────────────────────────────────────────────────────────────────────

# 运行时通过顶层 EMSuite 模块调用，避免循环引用
function _call_mie_pec(radius, freq, theta_vec)
    f = getfield(_emsuite_module(), :calculate_mie_rcs_pec_sphere)
    return f(radius, freq, theta_vec)
end

function _call_mie_dielectric(radius, freq, theta_vec, eps_r, mu_r)
    f = getfield(_emsuite_module(), :calculate_mie_rcs_dielectric_sphere)
    return f(radius, freq, theta_vec, eps_r, mu_r)
end

function _call_mie_dielectric_fullpol(radius, freq, theta_vec, eps_r, mu_r)
    return _call_mie_dielectric(radius, freq, theta_vec, eps_r, mu_r)
end

function _call_mie_pec_fullpol(radius, freq, theta_vec)
    f = getfield(_emsuite_module(), :calculate_mie_rcs_pec_sphere_fullpol)
    return f(radius, freq, theta_vec)
end

function _emsuite_module()
    return parentmodule(parentmodule(@__MODULE__))
end

function _expand_phi_obs(phi_obs::Real, n::Integer)
    return fill(float(phi_obs), n)
end

function _expand_phi_obs(phi_obs::AbstractVector, n::Integer)
    length(phi_obs) == n || throw(ArgumentError("phi_obs length must match theta_obs length"))
    return collect(float.(phi_obs))
end

function _bistatic_angles(theta_obs, phi_obs, theta_inc, phi_inc, polarization)
    k_hat = _spherical_unit_vector(theta_inc, phi_inc)

    pol_vec = collect(float.(polarization))
    pol_transverse = pol_vec .- dot(pol_vec, k_hat) .* k_hat
    pol_norm = norm(pol_transverse)
    pol_norm > sqrt(eps(Float64)) ||
        throw(ArgumentError("polarization must have a non-zero component transverse to the propagation direction"))

    e_parallel = pol_transverse ./ pol_norm
    e_perp = cross(k_hat, e_parallel)

    theta_scat = Vector{Float64}(undef, length(theta_obs))
    phi_local = Vector{Float64}(undef, length(theta_obs))
    @inbounds for i in eachindex(theta_obs)
        r_hat = _spherical_unit_vector(theta_obs[i], phi_obs[i])
        cosθ = clamp(dot(r_hat, k_hat), -1.0, 1.0)
        theta_scat[i] = acos(cosθ)

        transverse = r_hat .- cosθ .* k_hat
        sinθ = norm(transverse)
        if sinθ <= 1e-12
            phi_local[i] = 0.0
        else
            transverse ./= sinθ
            phi_local[i] = atan(dot(transverse, e_perp), dot(transverse, e_parallel))
        end
    end
    return theta_scat, phi_local
end

function _spherical_unit_vector(theta::Real, phi::Real)
    sθ, cθ = sincos(theta)
    sφ, cφ = sincos(phi)
    return [sθ * cφ, sθ * sφ, cθ]
end

# ─────────────────────────────────────────────────────────────────────────────
# 半波偶极子解析输入阻抗
# ─────────────────────────────────────────────────────────────────────────────

"""
    dipole_halfwave_Zin_analytic() -> ComplexF64

返回精确半波偶极子（L=λ/2）的输入阻抗解析值。

来源：Balanis《Antenna Theory》Table 8.1
```
Z_in ≈ 73.1 + j42.5  Ω
```
注：此值对应无限细导线。有限截面天线（如仿真中的薄圆柱）会有偏差，
实际验收标准：|ΔZ_in| / |Z_analytic| < 5%。
"""
function dipole_halfwave_Zin_analytic()
    return ComplexF64(73.1, 42.5)
end

"""
    dipole_resonant_Zin_analytic(; L_over_lambda = 0.47) -> ComplexF64

返回谐振偶极子的近似解析输入阻抗。

对 L ≈ 0.47λ，天线接近谐振，Im(Z_in) ≈ 0，通常 Re(Z_in) ≈ 73 Ω。
精确值依赖导线粗细，此处返回理想值作为参考上限。
"""
function dipole_resonant_Zin_analytic(; L_over_lambda::Float64 = 0.47)
    # 经验公式近似（Stutzman & Thiele）：谐振时 Z_in ≈ 73 + j0 Ω
    return ComplexF64(73.0, 0.0)
end

# ─────────────────────────────────────────────────────────────────────────────
# 半波偶极子解析远场方向图
# ─────────────────────────────────────────────────────────────────────────────

"""
    dipole_halfwave_farfield_analytic(theta_rad_vec) -> Vector{Float64}

计算半波偶极子的 E 面远场方向函数（归一化）。

## 公式（Balanis §8.3）
```
f(θ) = cos(π/2 · cos θ) / sin θ
```
沿 z 轴取向，入射极化沿 x，θ=π/2 时取最大值（赤道方向）。

## 返回
归一化方向函数（最大值 = 1.0），单位线性。

## 说明
- θ=0（z 轴方向）时 f=0（零点），计算时对 sin θ → 0 做保护
- 与 RCS 对比用途：乘以 |I₀|² 可得到辐射强度
"""
function dipole_halfwave_farfield_analytic(theta_rad_vec::AbstractVector)
    f = similar(Vector{Float64}, length(theta_rad_vec))
    for (i, θ) in enumerate(theta_rad_vec)
        sinθ = abs(sin(θ))
        if sinθ < 1e-10
            f[i] = 0.0
        else
            f[i] = abs(cos((π / 2) * cos(θ)) / sinθ)
        end
    end
    # 归一化：最大值 = 1.0
    mx = maximum(f)
    mx > 0 && (f ./= mx)
    return f
end

end # module ReferenceData
