"""
Verify SCFIE coupling sign convention against Legacy kernel formula.

Legacy kernel: Fsv ∝ (ρm·ρn/6 - 1/k²) × G_no4π
This test computes the coupling using both:
  A) EMSuite's scfie_coupling_interaction (current code)  
  B) The same quadrature but with the Legacy kernel formula
And checks if they match, and which sign convention (term1+term2 vs term1-term2) is correct.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations.SCFIEModule: scfie_coupling_interaction
using LinearAlgebra
using StaticArrays
using Printf

# Load mesh
mesh_file = joinpath(@__DIR__, "..", "..", "MoM_Kernels", "meshfiles", "TriTetra.nas")
tri_mesh, tet_mesh = read_mixed_nas_mesh(mesh_file, scale=0.001)

# Create basis
rwg = RWGBasis(tri_mesh)
swg = SWGBasis(tet_mesh)

# Frequency
freq = 2e9
k = 2π * freq / 3e8
omega = 2π * freq
mu0 = 4π * 1e-7
eps0 = 8.854187817e-12

# SCFIE
permittivities = [complex(2.0) for _ in 1:num_elements(tet_mesh)]
scfie = SCFIE(freq, permittivities)

# Get geometry info
tris = get_triangles_info(tri_mesh, rwg)
tetras = get_tetrahedra_info(tet_mesh, swg, permittivities)

# Pick elements (try a pair that aren't too close)
tri_info = tris[1]
tet_info = tetras[1]

println("=" ^ 60)
println("Coupling sign verification")
println("=" ^ 60)
println("κ = $(tet_info.κ)")
κ = tet_info.κ

# EMSuite's result
Z_sv_emsuite, Z_vs_emsuite = scfie_coupling_interaction(scfie, tri_info, tet_info)

# Now compute using Legacy-style kernel (ρ·ρ/6 - 1/k²) with proper EMSuite normalization
# This matches EMSuite's quadrature but uses the combined kernel from Legacy
gq_s = scfie.gq_surf
gq_v = scfie.gq_vol
Nq_s = length(gq_s.weight)
Nq_v = length(gq_v.weight)

r_q_s = tri_info.vertices * gq_s.coordinate
r_q_v = tet_info.vertices * gq_v.coordinate

vol_factor = tri_info.area * tet_info.volume

# Legacy combined kernel: jωμ₀ × (ρ·ρ/6 − 1/k²) × G₄π
# where ρ includes signed l/(2A) and A/(3V) normalization
# So: jωμ₀ × [(l/(2A)×ρ_m)·(A/(3V)×ρ_n) − (l/A)×(A/V)/(jωε₀×jωμ₀)] × G × vol_factor × w
# Wait, let me be explicit. The COMBINED Legacy kernel is:
#   jωμ₀ × f·f × G − jωμ₀/k² × div·div × G
# = jωμ₀ × f·f × G − 1/(jωε₀) × div·div × G     [since jωμ₀/k² = 1/(jωε₀)]
# = c1 * f·f * G + c2 * div·div * G               [c1=jωμ₀, c2=-1/(jωε₀) = c2 if using positive convention]

# Method: compute f·f and div·div per quadrature point and combine

# c1_legacy = jωμ₀ (without κ for now)
c1_leg = im * omega * mu0
# c2_legacy = 1/(jωε₀) = -j/(ωε₀)
c2_leg = 1.0 / (im * omega * eps0)

# Legacy formula: Z = c1 * <f,f'>_G + c2 * <div,div'>_G  (both terms, same sign)
# This comes from the L operator: Z = jωμ₀ <f,f'>_G + 1/(jωε₀) <div,div'>_G
# Which gives Z = jωμ₀ <f,f'>_G - j/(ωε₀) <div,div'>_G (scalar potential term negative)

Z_sv_plus = zeros(ComplexF64, 3, 4)   # term1 + term2 (Legacy-like)
Z_sv_minus = zeros(ComplexF64, 3, 4)  # term1 - term2 (current EMSuite)
Z_vs_plus = zeros(ComplexF64, 4, 3)  
Z_vs_minus = zeros(ComplexF64, 4, 3)  

for j in 1:Nq_v
    wj = gq_v.weight[j]
    rj = r_q_v[:, j]
    for i in 1:Nq_s
        wi = gq_s.weight[i]
        ri = r_q_s[:, i]
        R = norm(ri - rj)
        if R < 1e-10
            G = zero(ComplexF64)
        else
            G = exp(-im * k * R) / (4π * R)
        end
        factor = wi * wj * vol_factor * G

        for m in 1:3
            div_s = (tri_info.edgel[m] / tri_info.area) * tri_info.bfsSign[m]
            const_s = (tri_info.edgel[m] / (2 * tri_info.area)) * tri_info.bfsSign[m]
            f_m = const_s * (ri - tri_info.vertices[:, m])
            d_m = div_s
            
            for n in 1:4
                div_v = (tet_info.facesArea[n] / tet_info.volume) * tet_info.bfsSign[n]
                const_v = (tet_info.facesArea[n] / (3 * tet_info.volume)) * tet_info.bfsSign[n]
                f_n = const_v * (rj - tet_info.vertices[:, n])
                d_n = div_v
                
                ff = dot(f_m, f_n)
                dd = d_m * d_n
                
                # Z_SV with κ
                t1_sv = im * omega * mu0 * κ * ff
                t2_sv = κ / (im * omega * eps0) * dd
                Z_sv_plus[m, n]  += (t1_sv + t2_sv) * factor
                Z_sv_minus[m, n] += (t1_sv - t2_sv) * factor
                
                # Z_VS (no κ, but which c1 sign?)
                # Option A: c1_vs = +jωμ₀ (Legacy convention)
                t1_vs_pos = im * omega * mu0 * ff
                t2_vs = 1.0 / (im * omega * eps0) * dd
                Z_vs_plus[n, m]  += (t1_vs_pos + t2_vs) * factor
                Z_vs_minus[n, m] += (t1_vs_pos - t2_vs) * factor
            end
        end
    end
end

println("\n--- Z_SV comparison ---")
println("EMSuite uses (term1 - term2) with c1_sv = +jωμ₀κ")
println()
for m in 1:3, n in 1:4
    e = Z_sv_emsuite[m, n]
    p = Z_sv_plus[m, n]
    mi = Z_sv_minus[m, n]
    @printf("[%d,%d] EMSuite: %+.4e%+.4ej  Plus: %+.4e%+.4ej  Minus: %+.4e%+.4ej\n",
            m, n, real(e), imag(e), real(p), imag(p), real(mi), imag(mi))
end

# Check which variant matches
err_plus_sv = norm(Z_sv_emsuite - Z_sv_plus) / norm(Z_sv_emsuite)
err_minus_sv = norm(Z_sv_emsuite - Z_sv_minus) / norm(Z_sv_emsuite)
@printf("\nZ_SV: err(term1+term2) = %.2e, err(term1-term2) = %.2e\n", err_plus_sv, err_minus_sv)
if err_minus_sv < 1e-10
    println("→ EMSuite matches (term1 - term2) ✓ [as coded]")
elseif err_plus_sv < 1e-10
    println("→ EMSuite matches (term1 + term2) ← BUG: should use + but uses -")
end

println("\n--- Z_VS comparison ---")
println("EMSuite uses (term1_vs + term2_vs) with c1_vs = -jωμ₀")
println()
for n in 1:4, m in 1:3
    e = Z_vs_emsuite[n, m]
    p = Z_vs_plus[n, m]
    mi = Z_vs_minus[n, m]
    @printf("[%d,%d] EMSuite: %+.4e%+.4ej  Plus(+c1): %+.4e%+.4ej  Minus(+c1): %+.4e%+.4ej\n",
            n, m, real(e), imag(e), real(p), imag(p), real(mi), imag(mi))
end

err_plus_vs = norm(Z_vs_emsuite - Z_vs_plus) / norm(Z_vs_emsuite)
err_minus_vs = norm(Z_vs_emsuite - Z_vs_minus) / norm(Z_vs_emsuite)

# Also check with -c1 (current code has c1_vs = -jωμ₀)
Z_vs_negc1_plus = zeros(ComplexF64, 4, 3)
Z_vs_negc1_minus = zeros(ComplexF64, 4, 3)
for j in 1:Nq_v
    wj = gq_v.weight[j]
    rj = r_q_v[:, j]
    for i in 1:Nq_s
        wi = gq_s.weight[i]
        ri = r_q_s[:, i]
        R = norm(ri - rj)
        R < 1e-10 && continue
        G = exp(-im * k * R) / (4π * R)
        factor = wi * wj * vol_factor * G
        for m in 1:3
            const_s = (tri_info.edgel[m] / (2 * tri_info.area)) * tri_info.bfsSign[m]
            f_m = const_s * (ri - tri_info.vertices[:, m])
            d_m = (tri_info.edgel[m] / tri_info.area) * tri_info.bfsSign[m]
            for n in 1:4
                const_v = (tet_info.facesArea[n] / (3 * tet_info.volume)) * tet_info.bfsSign[n]
                f_n = const_v * (rj - tet_info.vertices[:, n])
                d_n = (tet_info.facesArea[n] / tet_info.volume) * tet_info.bfsSign[n]
                
                t1 = -im * omega * mu0 * dot(f_n, f_m)
                t2 = 1.0 / (im * omega * eps0) * d_n * d_m
                Z_vs_negc1_plus[n, m]  += (t1 + t2) * factor
                Z_vs_negc1_minus[n, m] += (t1 - t2) * factor
            end
        end
    end
end

err_negc1_plus_vs = norm(Z_vs_emsuite - Z_vs_negc1_plus) / norm(Z_vs_emsuite)
err_negc1_minus_vs = norm(Z_vs_emsuite - Z_vs_negc1_minus) / norm(Z_vs_emsuite)

@printf("\nZ_VS with c1=+jωμ₀: err(+) = %.2e, err(-) = %.2e\n", err_plus_vs, err_minus_vs)
@printf("Z_VS with c1=-jωμ₀: err(+) = %.2e, err(-) = %.2e\n", err_negc1_plus_vs, err_negc1_minus_vs)

# Determine which option the EMSuite code currently matches
options = [
    ("c1=+jωμ₀, term1+term2", err_plus_vs),
    ("c1=+jωμ₀, term1-term2", err_minus_vs),
    ("c1=-jωμ₀, term1+term2", err_negc1_plus_vs),
    ("c1=-jωμ₀, term1-term2", err_negc1_minus_vs)
]
best = argmin([o[2] for o in options])
println("\n→ EMSuite Z_VS matches: $(options[best][1]) (err = $(options[best][2]))")

println("\n" * "=" ^ 60)
println("EXPECTED (Legacy): both Z_SV and Z_VS use (term1 + term2) with c1 = +jωμ₀")
println("Legacy sign convention: Z_SV = κ × L(f_surf, f_vol), Z_VS = L(f_vol, f_surf)")
println("where L = jωμ₀ × <f,f'>_G + 1/(jωε₀) × <div,div'>_G")
println("=" ^ 60)

# Final check: does Z_SV/κ ≈ Z_VS^T (reciprocity)?
if !iszero(κ)
    Z_sv_normalized = Z_sv_emsuite ./ κ
    Z_vs_transposed = transpose(Z_vs_emsuite)
    recip_err = norm(Z_sv_normalized - Z_vs_transposed) / norm(Z_sv_normalized)
    @printf("\nReciprocity check: |Z_SV/κ - Z_VS^T| / |Z_SV/κ| = %.2e\n", recip_err)
    if recip_err < 1e-10
        println("✓ Reciprocity holds")
    else
        println("✗ Reciprocity violated! (indicates sign bug)")
    end
end
