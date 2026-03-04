"""
plot_rcs_sphere.jl — Sphere Bistatic RCS: Full-sphere Comparison

Validation against Mie series on a complete (θ, φ) sphere grid, not just
E/H-plane cuts. For incident +z, x-polarized wave the exact Mie prediction is:

    σ_θθ(θ, φ) = |S₂(θ)|² (4π/k²) · cos²φ
    σ_φφ(θ, φ) = |S₁(θ)|² (4π/k²) · sin²φ
    σ_tot(θ, φ) = σ_θθ + σ_φφ

Cases:
  1. PEC sphere  — EFIE + RWG     vs Mie (calculate_mie_rcs_pec_sphere_fullpol)
  2. Dielectric  — PMCHW + RWG    vs Mie (calculate_mie_rcs_dielectric_sphere)

Output plots (2×2 layout):
  [1,1] PEC   E/H-plane RCS curves + Legacy reference
  [1,2] PEC   full-sphere signed error heatmap (MoM − Mie) [dB]
  [2,1] PMCHW E/H-plane RCS curves
  [2,2] PMCHW full-sphere signed error heatmap [dB]

Usage:
  julia --project scripts/plot_rcs_sphere.jl
"""

using EMSuite
using Plots
using LinearAlgebra, Statistics, Printf, DelimitedFiles

gr()   # GR backend: no GUI, writes PNG

# ─────────────────────────────────────────────────────────────────
# Physical parameters
# ─────────────────────────────────────────────────────────────────
radius = 0.5    # m  (ka ≈ π at 300 MHz)
freq   = 300e6  # Hz

# ─────────────────────────────────────────────────────────────────
# Full-sphere observation grid (5-degree resolution)
# theta: 0→π  (Nth points)   phi: 0→2π (Nph points)
# ─────────────────────────────────────────────────────────────────
Nth        = 37    # 0:5:180  (37 points)
Nph        = 73    # 0:5:360  (73 points)
theta_grid = collect(range(0.0, π,  Nth))
phi_grid   = collect(range(0.0, 2π, Nph))
theta_d    = rad2deg.(theta_grid)
phi_d      = rad2deg.(phi_grid)

# sin(θ)-weighting: accounts for solid-angle dΩ ∝ sin(θ)dθdφ
# poles (sin→0) are set to zero to avoid pole artifacts
sin_w      = sin.(theta_grid)
sin_w[1]   = 0.0
sin_w[end] = 0.0

# Helper: solid-angle-weighted RMS error over full sphere
function full_sphere_rms(diff_db::Matrix{Float64})
    W = sin_w .* ones(1, size(diff_db, 2))   # Nth × Nph
    return sqrt(sum(W .* diff_db .^ 2) / (sum(sin_w) * size(diff_db, 2) + eps()))
end

# Convert linear RCS to dB, flooring at -40 dBsm to avoid -Inf
to_db(x) = 10.0 .* log10.(max.(x, 1e-4))

# ─────────────────────────────────────────────────────────────────
# Case 1 — PEC sphere (EFIE)
# ─────────────────────────────────────────────────────────────────
println("=" ^ 60)
println("Case 1: PEC Sphere (EFIE)")
println("=" ^ 60)

mesh_pec  = generate_sphere_mesh(radius, 16, 32)
basis_pec = RWGBasis(mesh_pec)
efie      = EFIE(freq)
Z_pec     = assemble_impedance_matrix(efie, basis_pec)
pw_pec    = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])   # +z propagation, x-polarized
V_pec     = excitation_vector(efie, pw_pec, basis_pec)
I_pec     = Z_pec \ V_pec
println("  Solved.  N = $(num_basis(basis_pec))")

set_frequency!(freq)
println("  Computing far field on $(Nth)×$(Nph) sphere grid ...")
rcs_pec_full, _, _ = radarCrossSection(theta_grid, phi_grid, I_pec, basis_pec)
# rcs_pec_full[comp, ith, iph]   comp=1 → θθ,  comp=2 → φφ
mom_tot_pec = rcs_pec_full[1, :, :] .+ rcs_pec_full[2, :, :]   # Nth × Nph

# Mie full-sphere: σ_tot(θ,φ) = |S₂(θ)|² cos²φ + |S₁(θ)|² sin²φ
mie_S2_pec, mie_S1_pec = calculate_mie_rcs_pec_sphere_fullpol(radius, freq, theta_grid)
cos2ph     = cos.(phi_grid)' .^ 2    # 1 × Nph (broadcast over theta)
sin2ph     = sin.(phi_grid)' .^ 2
mie_tot_pec = mie_S2_pec .* cos2ph .+ mie_S1_pec .* sin2ph   # Nth × Nph

mom_tot_pec_db = to_db(mom_tot_pec)
mie_tot_pec_db = to_db(mie_tot_pec)
err_pec_db     = mom_tot_pec_db .- mie_tot_pec_db   # signed: + means MoM > Mie

rms_pec    = full_sphere_rms(err_pec_db)
maxabs_pec = maximum(abs.(err_pec_db))
@printf "  Full-sphere weighted RMS = %.2f dB   max |err| = %.2f dB\n" rms_pec maxabs_pec

# E-plane (phi=0, index 1)
E_mom_pec = mom_tot_pec_db[:, 1];     E_mie_pec = mie_tot_pec_db[:, 1]
@printf "  E-plane RMS (phi=0)  = %.2f dB\n" sqrt(mean((E_mom_pec .- E_mie_pec).^2))

# H-plane (phi=π/2, nearest index)
H_idx     = argmin(abs.(phi_grid .- π/2))
H_mom_pec = mom_tot_pec_db[:, H_idx]; H_mie_pec = mie_tot_pec_db[:, H_idx]
@printf "  H-plane RMS (phi=90) = %.2f dB\n" sqrt(mean((H_mom_pec .- H_mie_pec).^2))

# Legacy reference (allinone_rcs.txt: theta [deg], RCS [dBsm])
legacy_rcs_path = joinpath(@__DIR__, "..", "..", "..", "LegacyBenchmark", "allinone_rcs.txt")
legacy_theta_d = nothing; legacy_rcs_db = nothing
if isfile(legacy_rcs_path)
    raw = readdlm(legacy_rcs_path)
    legacy_theta_d = raw[:, 1]
    legacy_rcs_db  = raw[:, 2]
    L = min(Nth, size(raw, 1))
    @printf "  Legacy E-plane RMS   = %.2f dB\n" sqrt(mean((E_mom_pec[1:L] .- legacy_rcs_db[1:L]).^2))
else
    println("  [INFO] Legacy allinone_rcs.txt not found, skipping overlay")
end

# ─────────────────────────────────────────────────────────────────
# Case 2 — Dielectric sphere (PMCHW)
# ─────────────────────────────────────────────────────────────────
println()
println("=" ^ 60)
println("Case 2: Dielectric Sphere (PMCHW)  eps_r=4")
println("=" ^ 60)

eps_r = 4.0; mu_r = 1.0
mesh_die  = generate_sphere_mesh(radius, 16, 32)
basis_die = RWGBasis(mesh_die)
pmchw     = PMCHW(freq, eps_r, mu_r)
Z_die     = assemble_impedance_matrix(pmchw, basis_die)
pw_die    = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
V_die     = excitation_vector(pmchw, pw_die, basis_die)
I_die     = Z_die \ V_die
println("  Solved.  N = $(num_basis(basis_die))")

println("  Computing far field on $(Nth)×$(Nph) sphere grid ...")
rcs_die_full, _, _ = radarCrossSection(theta_grid, phi_grid, I_die, basis_die,
                                        pmchw.k0, pmchw.eta0)
mom_tot_die = rcs_die_full[1, :, :] .+ rcs_die_full[2, :, :]   # Nth × Nph

# Mie full-sphere for dielectric: rcs_E → |S₂|², rcs_H → |S₁|²
mie_S2_die, mie_S1_die, _ = calculate_mie_rcs_dielectric_sphere(radius, freq,
                                                                   theta_grid, eps_r, mu_r)
mie_tot_die = mie_S2_die .* cos2ph .+ mie_S1_die .* sin2ph   # Nth × Nph

mom_tot_die_db = to_db(mom_tot_die)
mie_tot_die_db = to_db(mie_tot_die)
err_die_db     = mom_tot_die_db .- mie_tot_die_db

rms_die    = full_sphere_rms(err_die_db)
maxabs_die = maximum(abs.(err_die_db))
@printf "  Full-sphere weighted RMS = %.2f dB   max |err| = %.2f dB\n" rms_die maxabs_die

E_mom_die = mom_tot_die_db[:, 1];      E_mie_die = mie_tot_die_db[:, 1]
H_mom_die = mom_tot_die_db[:, H_idx];  H_mie_die = mie_tot_die_db[:, H_idx]
@printf "  E-plane RMS = %.2f dB   H-plane RMS = %.2f dB\n" sqrt(mean((E_mom_die.-E_mie_die).^2)) sqrt(mean((H_mom_die.-H_mie_die).^2))

# ─────────────────────────────────────────────────────────────────
# Plots — 2×2 layout
# ─────────────────────────────────────────────────────────────────
println("\nGenerating plots ...")

# [1,1]: PEC E/H-plane curves
p11 = plot(theta_d, E_mie_pec;
    label = "Mie E-plane", lw = 2, lc = :black, ls = :dash,
    xlabel = "Theta (deg)", ylabel = "RCS (dBsm)",
    title  = "PEC Sphere RCS (EFIE vs Mie)\nr=0.5m, f=300MHz",
    legend = :bottomright, xticks = 0:30:180, ylims = (-30, 10))
plot!(p11, theta_d, H_mie_pec;
    label = "Mie H-plane", lw = 2, lc = :gray, ls = :dash)
plot!(p11, theta_d, E_mom_pec;
    label = "EFIE E-plane (N=$(num_basis(basis_pec)))", lw = 1.5, lc = :steelblue)
plot!(p11, theta_d, H_mom_pec;
    label = "EFIE H-plane", lw = 1.5, lc = :royalblue, ls = :dot)
if !isnothing(legacy_theta_d)
    L = min(Nth, length(legacy_theta_d))
    plot!(p11, legacy_theta_d[1:L], legacy_rcs_db[1:L];
        label = "Legacy EFIE", lw = 1.5, lc = :darkorange, ls = :dot)
end
annotate!(p11, 90, -27,
    text("Full-sphere RMS = $(@sprintf("%.2f", rms_pec)) dB", 9, :gray))

# [1,2]: PEC error heatmap (signed, Nth × Nph)
err_clim_pec = max(1.0, ceil(maxabs_pec))
p12 = heatmap(phi_d, theta_d, err_pec_db;
    xlabel = "Phi (deg)", ylabel = "Theta (deg)",
    title  = "PEC Error (MoM−Mie) [dB]  RMS=$(round(rms_pec, digits=2)) dB",
    color  = :RdBu_11, clims = (-err_clim_pec, err_clim_pec),
    xticks = 0:90:360, yticks = 0:30:180)

# [2,1]: PMCHW E/H-plane curves
p21 = plot(theta_d, E_mie_die;
    label = "Mie E-plane", lw = 2, lc = :black, ls = :dash,
    xlabel = "Theta (deg)", ylabel = "RCS (dBsm)",
    title  = "Dielectric Sphere RCS (PMCHW vs Mie)\nr=0.5m, f=300MHz, eps_r=4",
    legend = :bottomright, xticks = 0:30:180, ylims = (-30, 10))
plot!(p21, theta_d, H_mie_die;
    label = "Mie H-plane", lw = 2, lc = :gray, ls = :dash)
plot!(p21, theta_d, E_mom_die;
    label = "PMCHW E-plane (N=$(num_basis(basis_die)))", lw = 1.5, lc = :tomato)
plot!(p21, theta_d, H_mom_die;
    label = "PMCHW H-plane", lw = 1.5, lc = :firebrick, ls = :dot)
annotate!(p21, 90, -27,
    text("Full-sphere RMS = $(@sprintf("%.2f", rms_die)) dB", 9, :gray))

# [2,2]: PMCHW error heatmap
err_clim_die = max(1.0, ceil(maxabs_die))
p22 = heatmap(phi_d, theta_d, err_die_db;
    xlabel = "Phi (deg)", ylabel = "Theta (deg)",
    title  = "Dielectric Error (MoM−Mie) [dB]  RMS=$(round(rms_die, digits=2)) dB",
    color  = :RdBu_11, clims = (-err_clim_die, err_clim_die),
    xticks = 0:90:360, yticks = 0:30:180)

# ─────────────────────────────────────────────────────────────────
# Combine and Save
# ─────────────────────────────────────────────────────────────────
combined = plot(p11, p12, p21, p22;
    layout = (2, 2), size = (1200, 900), dpi = 110,
    left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)

out_path = joinpath(dirname(@__DIR__), "docs", "images", "rcs_sphere_comparison.png")
savefig(combined, out_path)
println("\nImage saved: $out_path")
println()
println("Summary:")
@printf "  PEC  EFIE  full-sphere RMS = %.2f dB   max = %.2f dB\n" rms_pec  maxabs_pec
@printf "  PMCHW      full-sphere RMS = %.2f dB   max = %.2f dB\n" rms_die  maxabs_die
