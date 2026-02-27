using EMSuite
using JLD2
using LinearAlgebra
using Test
using StaticArrays
using Printf

# Helper to compute area
function compute_triangle_area(mesh, tri_idx)
    v_ids = mesh.triangles[:, tri_idx]
    v1 = mesh.node[:, v_ids[1]]
    v2 = mesh.node[:, v_ids[2]]
    v3 = mesh.node[:, v_ids[3]]
    return 0.5 * norm(cross(v2 - v1, v3 - v1))
end

# Load Golden Data
golden_file = joinpath(@__DIR__, "..", "legacy_golden_data.jld2")
if !isfile(golden_file)
    error("Golden data not found. Run generate_legacy_golden_data.jl first.")
end
@load golden_file Z bf_data tri_data frequency

println("Loaded Golden Data: $(length(bf_data)) BFs, $(length(tri_data)) Triangles.")

# Load Mesh using EMSuite
mesh_file = joinpath(@__DIR__, "..", "temp_sphere.nas")
mesh = EMSuite.Geometry.read_nas_mesh(mesh_file)

# Generate Basis Functions using EMSuite
println("Generating EMSuite Basis Functions...")
rwg_basis = EMSuite.BasisFunctions.RWGBasis(mesh)
ems_bfs = rwg_basis.functions

println("EMSuite Generated: $(length(ems_bfs)) BFs.")

# Verification Loop
println("\n--- Verifying Basis Functions ---")

@testset "Basis Function Verification" begin
    @test length(ems_bfs) == length(bf_data)
    
    mismatches = 0
    max_len_diff = 0.0
    max_center_diff = 0.0
    
    # We assume the order is the same because we tried to match the sorting logic.
    # If there are massive mismatches, we might need to build a map.
    
    for i in 1:length(ems_bfs)
        ems_bf = ems_bfs[i]
        leg_bf = bf_data[i] # Dict
        
        # 1. Check Edge Length
        diff_len = abs(ems_bf.edge_length - leg_bf["edge_length"])
        max_len_diff = max(max_len_diff, diff_len)
        if diff_len > 1e-10
            println("BF $i Length Mismatch: EMS=$(ems_bf.edge_length), LEG=$(leg_bf["edge_length"])")
            mismatches += 1
        end
        
        # 2. Check Support Triangles
        # Legacy: inGeo (Vector of 2)
        # EMSuite: support (SVector of 2)
        ems_tris = sort([ems_bf.support[1], ems_bf.support[2]])
        leg_tris = sort([leg_bf["triangles"][1], leg_bf["triangles"][2]])
        
        if ems_tris != leg_tris
            println("BF $i Support Mismatch: EMS=$ems_tris, LEG=$leg_tris")
            mismatches += 1
        end
        
        # 3. Check Center
        diff_center = norm(ems_bf.center - leg_bf["center"])
        max_center_diff = max(max_center_diff, diff_center)
        if diff_center > 1e-10
            println("BF $i Center Mismatch: EMS=$(ems_bf.center), LEG=$(leg_bf["center"])")
            mismatches += 1
        end
    end
    
    println("Total BF Mismatches: $mismatches")
    println("Max Length Diff: $max_len_diff")
    println("Max Center Diff: $max_center_diff")
    
    @test mismatches == 0
end

println("\n--- Verifying Triangle Geometry ---")
@testset "Triangle Geometry Verification" begin
    @test EMSuite.Geometry.num_elements(mesh) == length(tri_data)
    
    max_area_diff = 0.0
    
    for i in 1:length(tri_data)
        leg_tri = tri_data[i]
        ems_area = compute_triangle_area(mesh, i)
        leg_area = leg_tri["area"]
        
        diff_area = abs(ems_area - leg_area)
        max_area_diff = max(max_area_diff, diff_area)
        
        if diff_area > 1e-10
            println("Tri $i Area Mismatch: EMS=$ems_area, LEG=$leg_area")
        end
    end
    
    println("Max Area Diff: $max_area_diff")
    @test max_area_diff < 1e-10
end
