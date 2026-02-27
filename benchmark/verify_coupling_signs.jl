"""
Verify SCFIE coupling (Z_SV, Z_VS) element-level signs against Legacy formula.

Legacy formula (from MoM_Kernels/src/ZmatAndVvec/EFIE/EFIEVSIERWGSWG.jl):
    Fsv = lman * Σ((ρm·ρn/6 - divk²) * gw)
    Z_SV = κ * JKη_0div4π * Fsv
    Z_VS = JKη_0div4π * Fsv

Where:
    JKη_0div4π = jωμ₀/(4π)
    divk² = 1/k²
    lman = signed_lm * signed_An
    gw = exp(-jkR)/R * w_i * w_j
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
if !isfile(mesh_file)
    error("Mesh file not found: $mesh_file")
end

mesh = read_mixed_nas_mesh(mesh_file, scale=0.001)  # mm → m
tri_mesh, tet_mesh = mesh

# Create basis functions
rwg = RWGBasis(tri_mesh)
swg = SWGBasis(tet_mesh)

# Frequency
freq = 2e9  # 2 GHz
k = 2π * freq / 3e8
omega = 2π * freq
mu0 = 4π * 1e-7
eps0 = 8.854187817e-12
eta = sqrt(mu0 / eps0)

# Create SCFIE operator (needs permittivities)
# TriTetra mesh has dielectric tetra with εr = 2.0 typically
permittivities = [complex(2.0) for _ in 1:num_elements(tet_mesh)]
scfie = SCFIE(freq, permittivities)

# Pick the first triangle and first tetrahedron
tris = get_triangles_info(tri_mesh, rwg)
tetras = get_tetrahedra_info(tet_mesh, swg, permittivities)

tri_info = tris[1]
tet_info = tetras[1]

# Compute coupling using EMSuite's function
Z_sv_emsuite, Z_vs_emsuite = scfie_coupling_interaction(scfie, tri_info, tet_info)

println("=" ^ 60)
println("Element-level coupling sign verification")
println("=" ^ 60)
println("Tri idx = 1, Tet idx = 1")
println("κ = $(tet_info.κ)")
println("tri.area = $(tri_info.area)")
println("tet.volume = $(tet_info.volume)")
println("tri.edgel = $(tri_info.edgel)")
println("tet.facesArea = $(tet_info.facesArea)")
println("tri.bfsSign = $(tri_info.bfsSign)")
println("tet.bfsSign = $(tet_info.bfsSign)")
println("vol_factor (area*vol) = $(tri_info.area * tet_info.volume)")
println()

# Compute coupling using Legacy's formula
# Legacy Green's function: G_legacy = exp(-jkR) / R  (no 1/(4π))
# Legacy factor: JKη_0div4π = jωμ₀/(4π)
# Legacy formula: Fsv = lman * Σ((ρm·ρn/6 - 1/k²) * exp(-jkR)/R * w_i * w_j)

JKeta0div4pi = im * omega * mu0 / (4π)
divk2 = 1.0 / k^2
kappa = tet_info.κ

# Quadrature 
gq_s = scfie.gq_surf
gq_v = scfie.gq_vol
Nq_s = length(gq_s.weight)
Nq_v = length(gq_v.weight)

# Quadrature points in physical space
r_q_s = tri_info.vertices * gq_s.coordinate  # 3 × Nq_s
r_q_v = tet_info.vertices * gq_v.coordinate  # 3 × Nq_v

# Compute Z_sv_legacy and Z_vs_legacy
Z_sv_legacy = zeros(ComplexF64, 3, 4)
Z_vs_legacy = zeros(ComplexF64, 4, 3)

for m in 1:3  # RWG basis on triangle
    # Signed edge length (Legacy absorbs sign into lm)
    signed_lm = tri_info.edgel[m] * tri_info.bfsSign[m]
    v_free_m = tri_info.vertices[:, m]
    
    for n in 1:4  # SWG basis on tetrahedron
        # Signed face area (Legacy absorbs sign into An)
        signed_An = tet_info.facesArea[n] * tet_info.bfsSign[n]
        v_free_n = tet_info.vertices[:, n]
        
        lman = signed_lm * signed_An
        
        # Double quadrature
        Fsv = zero(ComplexF64)
        for j in 1:Nq_v
            wj = gq_v.weight[j]
            rj = r_q_v[:, j]
            rho_n = rj - v_free_n
            
            for i in 1:Nq_s
                wi = gq_s.weight[i]
                ri = r_q_s[:, i]
                rho_m = ri - v_free_m
                
                R = norm(ri - rj)
                if R < 1e-10
                    G_legacy = zero(ComplexF64)
                else
                    G_legacy = exp(-im * k * R) / R
                end
                
                Fsv += (dot(rho_m, rho_n) / 6 - divk2) * G_legacy * wi * wj
            end
        end
        
        Fsv *= lman
        
        # Legacy Z_SV: κ * JKη_0div4π * Fsv
        Z_sv_legacy[m, n] = kappa * JKeta0div4pi * Fsv
        
        # Legacy Z_VS: JKη_0div4π * Fsv (no κ)
        Z_vs_legacy[n, m] = JKeta0div4pi * Fsv
    end
end

# Note: Legacy includes vol_factor (area * volume) implicitly through quadrature:
# Physical integration = area * volume * Σ(w_i * w_j * f(r_i, r_j))
# But our lman and ρ·ρ/6 implicitly cancel these factors.
# Need to add the area * volume factor
Z_sv_legacy .*= tri_info.area * tet_info.volume
Z_vs_legacy .*= tri_info.area * tet_info.volume

println("Z_SV comparison (3×4 matrix):")
println("-" ^ 60)
max_err_sv = 0.0
for m in 1:3, n in 1:4
    e = Z_sv_emsuite[m, n]
    l = Z_sv_legacy[m, n]
    ratio = iszero(l) ? NaN : e / l
    err = abs(e - l) / max(abs(e), abs(l), 1e-30)
    global max_err_sv = max(max_err_sv, err)
    @printf("  [%d,%d] EMSuite: %+.6e%+.6ej  Legacy: %+.6e%+.6ej  ratio: %.4f\n", 
            m, n, real(e), imag(e), real(l), imag(l), abs(ratio))
end
println()
@printf("Max relative error Z_SV: %.6e\n\n", max_err_sv)

println("Z_VS comparison (4×3 matrix):")
println("-" ^ 60)
max_err_vs = 0.0
for n in 1:4, m in 1:3
    e = Z_vs_emsuite[n, m]
    l = Z_vs_legacy[n, m]
    ratio = iszero(l) ? NaN : e / l
    err = abs(e - l) / max(abs(e), abs(l), 1e-30)
    global max_err_vs = max(max_err_vs, err)
    @printf("  [%d,%d] EMSuite: %+.6e%+.6ej  Legacy: %+.6e%+.6ej  ratio: %.4f\n", 
            n, m, real(e), imag(e), real(l), imag(l), abs(ratio))
end
println()
@printf("Max relative error Z_VS: %.6e\n\n", max_err_vs)

# Summary
println("=" ^ 60)
println("SUMMARY")
println("=" ^ 60)
if max_err_sv < 1e-10 && max_err_vs < 1e-10
    println("✓ Both Z_SV and Z_VS match Legacy formula!")
elseif max_err_sv < 1e-10
    println("✓ Z_SV matches Legacy")
    println("✗ Z_VS does NOT match Legacy (err = $max_err_vs)")
elseif max_err_vs < 1e-10
    println("✗ Z_SV does NOT match Legacy (err = $max_err_sv)")
    println("✓ Z_VS matches Legacy")
else
    println("✗ Z_SV does NOT match Legacy (err = $max_err_sv)")
    println("✗ Z_VS does NOT match Legacy (err = $max_err_vs)")
end
