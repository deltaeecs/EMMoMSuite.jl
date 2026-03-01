using EMSuite

# Check Tri.nas
m1 = read_nas_mesh(joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "Tri.nas"), scale=1.0)
println("Tri.nas: $(typeof(m1))")
if m1 isa EMSuite.Geometry.TriangleMesh
    b1 = RWGBasis(m1)
    println("  N = $(num_basis(b1))")
end

# Check TriEFIE.nas
m2 = read_nas_mesh(joinpath(@__DIR__, "..", "..", "MoM_AllinOne", "meshfiles", "TriEFIE.nas"), scale=1.0)
println("TriEFIE.nas: $(typeof(m2))")
if m2 isa EMSuite.Geometry.TriangleMesh
    b2 = RWGBasis(m2)
    println("  N = $(num_basis(b2))")
end

# Check plate_benchmark.nas
m3 = read_nas_mesh(joinpath(@__DIR__, "plate_benchmark.nas"), scale=1.0)
println("plate_benchmark.nas: $(typeof(m3))")
if m3 isa EMSuite.Geometry.TriangleMesh
    b3 = RWGBasis(m3)
    println("  N = $(num_basis(b3))")
end
