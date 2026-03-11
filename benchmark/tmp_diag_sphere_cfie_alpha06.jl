using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using Printf

const ROOT_DIR = joinpath(@__DIR__, "..")
const MESH_DIR = joinpath(ROOT_DIR, "..", "MoM_AllinOne", "meshfiles")
const FEKO_DIR = joinpath(ROOT_DIR, "..", "MoM_AllinOne", "deps", "compare_feko")

freq = 6.0e8
alpha = 0.6
set_frequency!(freq)
mesh = read_nas_mesh(joinpath(MESH_DIR, "sphere_600MHz.nas"), scale = 1.0)
basis = RWGBasis(mesh)
source = PlaneWave(freq, π / 2, π, [0.0, 0.0, 1.0])
θs_obs = collect(range(-π, π, length = 721))
ϕs_obs = [0.0, π / 2]

theta_deg, phi_deg, _, rcs_ref = read_feko_rcs(joinpath(FEKO_DIR, "sphere_600MHzRCS.csv"))
cuts = split_phi_cuts(theta_deg, phi_deg, rcs_ref)
mie_phi0 = mie_pec_rcs_dBsm(joinpath(MESH_DIR, "sphere_600MHz.nas"), freq, θs_obs)

op = CFIE(freq, alpha)
Z = assemble_impedance_matrix(op, basis)
V = excitation_vector(op, source, basis)
I = Z \ V
_, _, rcs_dB = radarCrossSection(θs_obs, ϕs_obs, I, basis)

r_feko_0 = compute_rcs_accuracy(rcs_dB[:, 1], cuts[0.0].rcs_dBsm, cuts[0.0].theta, "phi0"; threshold = 2.0)
r_feko_90 = compute_rcs_accuracy(rcs_dB[:, 2], cuts[90.0].rcs_dBsm, cuts[90.0].theta, "phi90"; threshold = 2.0)
r_mie_0 = compute_rcs_accuracy(rcs_dB[:, 1], mie_phi0, rad2deg.(θs_obs), "mie"; threshold = 2.0)

@printf("alpha=%4.2f  FEKO(phi0,phi90)=(%.3f, %.3f)  Mie(phi0)=%.3f\n",
    alpha, r_feko_0.rmse_dB, r_feko_90.rmse_dB, r_mie_0.rmse_dB)