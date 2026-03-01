using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations

# Check TriTetra mesh size
mesh_file = joinpath(@__DIR__, "../MoM_Basics/meshfiles/TriTetra.nas")
surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001)
surf_basis = RWGBasis(surf_mesh); vol_basis = SWGBasis(vol_mesh)
n_surf = num_basis(surf_basis); n_vol = num_basis(vol_basis)
println("RWG: $n_surf, SWG: $n_vol, Total: $(n_surf+n_vol)")

# Check Tetra mesh size
mesh_file2 = joinpath(@__DIR__, "../MoM_AllinOne/meshfiles/Tetra.nas")
mesh2 = read_nas_mesh(mesh_file2)
if maximum(abs.(mesh2.node)) > 10.0; mesh2.node .*= 0.001; end
basis2 = SWGBasis(mesh2)
println("Tetra SWG: $(num_basis(basis2))")
