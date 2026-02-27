# Quick verification: near interaction ratio after fix
using EMSuite
using LinearAlgebra
using Printf

using EMSuite.IntegralEquations.EFIEModule: EFIE, calc_near_interaction!
using EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry: TriangleInfo, GaussQuadratureInfo, get_global_quad_points
using EMSuite.IntegralEquations.EFIEModule.Singularities: faceSingularityIgIvecg
using EMSuite.CoreModule: num_elements

const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne")
mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)
freq = 1e8
set_frequency!(freq)
basis = RWGBasis(mesh)
efie = EFIE(freq)

bf1 = basis.functions[1]
tri_plus = TriangleInfo(mesh, bf1.support[1])
tri_minus = TriangleInfo(mesh, bf1.support[2])

# EMSuite near (with fix applied)
Z_near_ems = zeros(ComplexF64, 3, 3)
calc_near_interaction!(Z_near_ems, efie, tri_plus, tri_minus)
Z_near_ems .*= efie.factor

# Legacy formula (manual)
Z_near_leg = zeros(ComplexF64, 3, 3)
gq_near = efie.gq_near
r_test = get_global_quad_points(tri_plus, gq_near)
w_test = gq_near.weight
C4divk2 = efie.C4divk2

for gi in 1:length(w_test)
    ri = r_test[gi]
    Ig, IvecSg = faceSingularityIgIvecg(
        ri, tri_minus.vertices, tri_minus.edgel, tri_minus.edgev̂,
        tri_minus.edgen̂, tri_minus.area, tri_minus.facen̂, efie.SSCg)
    
    for ni in 1:3
        rho_n = ri - tri_minus.vertices[:, ni]
        for mi in 1:3
            rho_m = ri - tri_plus.vertices[:, mi]
            Ztemp = ((dot(rho_m, rho_n) - C4divk2) * Ig - dot(rho_m, IvecSg)) * w_test[gi] / tri_minus.area
            Z_near_leg[mi, ni] += Ztemp
        end
    end
end
for ni in 1:3, mi in 1:3
    Z_near_leg[mi, ni] *= tri_plus.edgel[mi] * tri_minus.edgel[ni] * efie.factor
end

println("Near interaction ratio after fix:")
for m in 1:3, n in 1:3
    if abs(Z_near_leg[m,n]) > 1e-15
        ratio = Z_near_ems[m,n] / Z_near_leg[m,n]
        @printf("  [%d,%d] ratio = %.6f\n", m, n, abs(ratio))
    end
end
println("\nZ_near_ems[1,1] = $(Z_near_ems[1,1])")
println("Z_near_leg[1,1] = $(Z_near_leg[1,1])")
