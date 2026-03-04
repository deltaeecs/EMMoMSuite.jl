"""
plot_rcs_sphere.jl — Sphere Bistatic RCS Visualization

Comparison:
  - PEC sphere (EFIE + RWG) vs Mie analytical solution
  - Homogeneous dielectric sphere (PMCHW + RWG) vs Mie (E-plane & H-plane)

Outputs:
  - Plot 1: PEC sphere RCS (E-plane, theta = 0->180 deg)
  - Plot 2: Dielectric sphere PMCHW RCS (E-plane & H-plane)

Usage:
  julia --project scripts/plot_rcs_sphere.jl
"""

using EMSuite
using Plots
using LinearAlgebra, Statistics, Printf, DelimitedFiles

gr()   # GR backend: no GUI needed, outputs PNG

# ─────────────────────────────────────────────────────────────────
# Common Parameters
# ─────────────────────────────────────────────────────────────────
radius  = 0.5          # m  (electrical size ka ~= pi)
freq    = 300e6        # Hz
theta_v = collect(range(0.0, π, 181))   # 0->180 deg
theta_d = rad2deg.(theta_v)

# ─────────────────────────────────────────────────────────────────
# Plot 1: PEC Sphere — EFIE vs Mie
# ─────────────────────────────────────────────────────────────────
println("=" ^ 60)
println("Plot 1: PEC Sphere EFIE vs Mie")

mesh_pec   = generate_sphere_mesh(radius, 16, 32)
basis_pec  = RWGBasis(mesh_pec)
efie       = EFIE(freq)
Z_pec      = assemble_impedance_matrix(efie, basis_pec)
pw_pec     = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])   # +z propagation, x-polarized (matches Mie convention)
V_pec      = excitation_vector(efie, pw_pec, basis_pec)
I_pec      = Z_pec \ V_pec

set_frequency!(freq)
rcs_pec, _, rcs_pec_db = radarCrossSection(theta_v, [0.0], I_pec, basis_pec)
rcs_pec_Eplane = 10 .* log10.(rcs_pec[1, :, 1] .+ 1e-30)

mie_pec     = calculate_mie_rcs_pec_sphere(radius, freq, theta_v)
mie_pec_db  = 10 .* log10.(mie_pec .+ 1e-30)

err_rms = sqrt(mean((rcs_pec_Eplane .- mie_pec_db).^2))
@printf "  EFIE N=%d  RMS error=%.2f dB\n" num_basis(basis_pec) err_rms

# Load Legacy EFIE reference data (allinone_rcs.txt: theta [deg], RCS [dBsm])
legacy_rcs_path = joinpath(@__DIR__, "..", "..", "..", "LegacyBenchmark", "allinone_rcs.txt")
legacy_data = nothing
legacy_theta_d = nothing
legacy_rcs_db  = nothing
if isfile(legacy_rcs_path)
    raw = readdlm(legacy_rcs_path)
    legacy_theta_d = raw[:, 1]   # degrees
    legacy_rcs_db  = raw[:, 2]   # dBsm
    err_legacy = sqrt(mean((rcs_pec_Eplane .- legacy_rcs_db[1:length(rcs_pec_Eplane)]).^2))
    @printf "  vs Legacy (allinone_rcs)  RMS=%.2f dB\n" err_legacy
else
    println("  [INFO] Legacy reference file not found, skipping overlay")
end

p1 = plot(theta_d, mie_pec_db,
    label = "Mie (analytical)",
    lw = 2, lc = :black, ls = :dash,
    xlabel = "Observation Angle theta (deg)",
    ylabel = "RCS (dBsm)",
    title  = "PEC Sphere Bistatic RCS (EFIE vs Mie)\nradius=0.5 m, f=300 MHz",
    legend = :bottomright,
    xticks = 0:30:180,
    ylims  = (-30, 5),
)
plot!(p1, theta_d, rcs_pec_Eplane,
    label = "EFIE (N=$(num_basis(basis_pec)))",
    lw = 1.5, lc = :steelblue,
)
if !isnothing(legacy_theta_d)
    plot!(p1, legacy_theta_d, legacy_rcs_db,
        label = "Legacy EFIE (MoM_AllinOne)",
        lw = 1.5, lc = :darkorange, ls = :dot,
    )
end
annotate!(p1, 90, -27, text("RMS vs Mie = $(@sprintf(\"%.2f\", err_rms)) dB", 9, :gray))

# ─────────────────────────────────────────────────────────────────
# Plot 2: Dielectric Sphere — PMCHW vs Mie (E-plane & H-plane)
# ─────────────────────────────────────────────────────────────────
println("=" ^ 60)
println("Plot 2: Dielectric Sphere PMCHW vs Mie (E-plane & H-plane)")

eps_r = 4.0;  mu_r = 1.0
mesh_die  = generate_sphere_mesh(radius, 16, 32)
basis_die = RWGBasis(mesh_die)
pmchw     = PMCHW(freq, eps_r, mu_r)
Z_die     = assemble_impedance_matrix(pmchw, basis_die)
pw_die    = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
V_die     = excitation_vector(pmchw, pw_die, basis_die)
I_die     = Z_die \ V_die

rcs_Ep, _, _ = radarCrossSection(theta_v, [0.0],  I_die, basis_die,
                                  pmchw.k0, pmchw.eta0)
rcs_Hp, _, _ = radarCrossSection(theta_v, [π/2], I_die, basis_die,
                                  pmchw.k0, pmchw.eta0)
mom_Eplane = 10 .* log10.(rcs_Ep[1, :, 1] .+ 1e-30)
mom_Hplane = 10 .* log10.(rcs_Hp[2, :, 1] .+ 1e-30)

mie_E, mie_H, _ = calculate_mie_rcs_dielectric_sphere(radius, freq, theta_v, eps_r, mu_r)
mie_Edb = 10 .* log10.(mie_E .+ 1e-30)
mie_Hdb = 10 .* log10.(mie_H .+ 1e-30)

err_E = sqrt(mean((mom_Eplane .- mie_Edb).^2))
err_H = sqrt(mean((mom_Hplane .- mie_Hdb).^2))
@printf "  PMCHW N=%d  E-plane RMS=%.2f dB  H-plane RMS=%.2f dB\n" num_basis(basis_die) err_E err_H

p2 = plot(theta_d, mie_Edb,
    label = "Mie E-plane",
    lw = 2, lc = :black, ls = :dash,
    xlabel = "Observation Angle theta (deg)",
    ylabel = "RCS (dBsm)",
    title  = "Dielectric Sphere Bistatic RCS (PMCHW vs Mie)\nradius=0.5 m, f=300 MHz, eps_r=4",
    legend = :bottomright,
    xticks = 0:30:180,
    ylims  = (-30, 5),
)
plot!(p2, theta_d, mie_Hdb,
    label = "Mie H-plane",
    lw = 2, lc = :gray, ls = :dash,
)
plot!(p2, theta_d, mom_Eplane,
    label = "PMCHW E-plane (N=$(num_basis(basis_die)))",
    lw = 1.5, lc = :steelblue,
)
plot!(p2, theta_d, mom_Hplane,
    label = "PMCHW H-plane",
    lw = 1.5, lc = :tomato,
)
annotate!(p2, 90, -27, text("E-plane RMS=$(@sprintf(\"%.2f\",err_E)) dB  H-plane RMS=$(@sprintf(\"%.2f\",err_H)) dB", 9, :gray))

# ─────────────────────────────────────────────────────────────────
# Combine and Save
# ─────────────────────────────────────────────────────────────────
combined = plot(p1, p2, layout = (1, 2), size = (1100, 450), dpi = 120)
out_path = joinpath(dirname(@__DIR__), "docs", "images", "rcs_sphere_comparison.png")
savefig(combined, out_path)
println("\nImage saved: $out_path")
