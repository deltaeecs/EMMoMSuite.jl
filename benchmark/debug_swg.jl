using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions

nodes = [0.0 1.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0 1.0; 0.0 0.0 0.0 1.0 1.0]
elements = [1 2; 2 3; 3 4; 4 5]
tags = [1, 1]
mesh = TetrahedraMesh(2, nodes, elements, tags)
basis = SWGBasis(mesh)

println("num_basis = ", num_basis(basis))
for (i, f) in enumerate(basis.functions)
    println("BF $i: support=$(f.support), area=$(f.area), is_boundary=$(f.is_boundary), local_face=$(f.local_face_idx), signs=$(f.signs)")
end
