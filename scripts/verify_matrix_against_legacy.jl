using EMSuite
using JLD2
using LinearAlgebra
using Test
using StaticArrays
using Printf
using Statistics

# Load Golden Data
golden_file = joinpath(@__DIR__, "..", "legacy_golden_data.jld2")
if !isfile(golden_file)
    error("Golden data not found. Run generate_legacy_golden_data.jl first.")
end
@load golden_file Z bf_data tri_data frequency

println("Loaded Golden Data: Z size $(size(Z))")

# Load Mesh
mesh_file = joinpath(@__DIR__, "..", "temp_sphere.nas")
mesh = EMSuite.Geometry.read_nas_mesh(mesh_file)

# Generate Basis
rwg_basis = EMSuite.BasisFunctions.RWGBasis(mesh)

# Setup EFIE
efie = EMSuite.IntegralEquations.EFIE(frequency)

# Compute Z Matrix
println("Computing EMSuite Impedance Matrix...")
Z_ems = EMSuite.IntegralEquations.assemble_impedance_matrix(efie, rwg_basis)

# Compare
println("\n--- Verifying Impedance Matrix ---")
@testset "Matrix Verification" begin
    @test size(Z_ems) == size(Z)
    
    diff_Z = abs.(Z_ems - Z)
    max_diff = maximum(diff_Z)
    rel_diff = norm(Z_ems - Z) / norm(Z)
    
    println("Max Element Diff: $max_diff")
    println("Relative Norm Diff: $rel_diff")
    
    # Check diagonal elements specifically (Self-term)
    diag_diff = abs.(diag(Z_ems) - diag(Z))
    max_diag_diff = maximum(diag_diff)
    println("Max Diagonal Diff: $max_diag_diff")
    
    # Calculate Ratio
    # Avoid division by zero
    mask = abs.(Z) .> 1e-10 * maximum(abs.(Z))
    ratio = Z_ems[mask] ./ Z[mask]
    avg_ratio = mean(ratio)
    std_ratio = std(ratio)
    
    println("Average Ratio (EMS/LEG): $avg_ratio")
    println("Std Dev of Ratio: $std_ratio")
    
    # Check if ratio is close to 1 or some constant
    if abs(avg_ratio - 1.0) < 1e-2
        println("Result: MATCH (Ratio ~ 1.0)")
    elseif abs(abs(avg_ratio) - 4.0) < 1e-2
        println("Result: FACTOR 4 DETECTED")
    elseif abs(abs(avg_ratio) - 0.25) < 1e-2
        println("Result: FACTOR 1/4 DETECTED")
    else
        println("Result: MISMATCH (Unknown Factor)")
    end
    
    println("\n--- Detailed Analysis ---")
    println("Mean |Z|: $(mean(abs.(Z)))")
    
    # Diagonal
    diag_ems = diag(Z_ems)
    diag_leg = diag(Z)
    diag_ratio = diag_ems ./ diag_leg
    println("Diagonal Mean Ratio: $(mean(diag_ratio))")
    println("Diagonal Max Diff: $(maximum(abs.(diag_ems - diag_leg)))")
    
    # Off-Diagonal
    off_mask = ones(Bool, size(Z))
    for i in 1:size(Z,1); off_mask[i,i] = false; end
    off_ratio = Z_ems[off_mask] ./ Z[off_mask]
    println("Off-Diagonal Mean Ratio: $(mean(off_ratio))")
    
    # Check for sign flips
    sign_diff = angle.(Z_ems) .- angle.(Z)
    # Wrap to -pi, pi
    sign_diff = mod2pi.(sign_diff .+ π) .- π
    println("Mean Phase Diff: $(mean(abs.(sign_diff)))")
    
    # Debug Singular Values for Triangle 1
    println("\n--- Debug Singular Values (EMSuite) ---")
    # Need to construct TriangleInfo manually or use internal function
    # get_triangle_info is in Impedance module but not exported?
    # It is used in assemble_generic.
    # Let's reconstruct it manually for Tri 1.
    
    v_ids = mesh.triangles[:, 1]
    v1 = mesh.node[:, v_ids[1]]
    v2 = mesh.node[:, v_ids[2]]
    v3 = mesh.node[:, v_ids[3]]
    
    # Edges: v2-v3, v3-v1, v1-v2
    l1 = norm(v2 - v3)
    l2 = norm(v3 - v1)
    l3 = norm(v1 - v2)
    
    area = 0.5 * norm(cross(v2 - v1, v3 - v1))
    
    println("Triangle 1:")
    println("  Edges: $([l1, l2, l3])")
    println("  Area: $area")
    
    C4divk2 = efie.C4divk2
    println("  C4divk2: $C4divk2")
    
    sF1 = EMSuite.IntegralEquations.EFIEModule.Singularities.singularF1(l1, l2, l3)
    F1 = C4divk2 * sF1
    println("  sF1: $sF1")
    println("  F1: $F1")
    
    sF21 = EMSuite.IntegralEquations.EFIEModule.Singularities.singularF21(l1, l2, l3, area^2)
    println("  sF21: $sF21")
    
    diff = sF21 - F1
    println("  Diff (sF21 - F1): $diff")
    
    println("\n--- Z[1,1] Comparison ---")
    println("Legacy Z[1,1]: $(Z[1,1])")
    println("EMSuite Z[1,1]: $(Z_ems[1,1])")
    println("Ratio: $(Z_ems[1,1] / Z[1,1])")

    # --- Component Analysis for Z[1,1] ---
    println("\n--- Component Analysis for Z[1,1] ---")
    # Load Basis Data if not already loaded
    if !@isdefined(basis_data)
        # data is loaded in the outer scope, but we are inside a testset?
        # No, we are inside a loop over triangles? No, we are inside the testset block.
        # `data` variable from line 15 should be available if not shadowed.
        # But `data` was loaded at top level.
        # Let's assume `data` is available.
        # Wait, the error says `data` not defined.
        # Ah, `data` is defined at top level.
        # But `testset` introduces a new scope?
        # Yes.
        # We need to pass `data` or reload it.
        basis_data = load(golden_file)["bf_data"]
    end
    
    bf_idx = 1
    bf_info = basis_data[bf_idx]
    tri_ids = bf_info["triangles"]
    
    println("Basis Function 1 Triangles: ", tri_ids)
    
    # Get TriangleInfo objects
    # get_triangle_info requires basis as second argument?
    # No, that's for basis-specific info.
    # We want generic TriangleInfo.
    # It seems `get_triangle_info` is not what we want.
    # We want to construct TriangleInfo from mesh.
    # Let's look at `src/Geometry/MeshTypes.jl` or `src/Core/Types.jl`.
    # Actually, `RWGBasis` constructor creates `TriangleInfo` internally? No.
    # `EFIE` uses `TriangleInfo`.
    # `assemble_generic` creates it.
    # Let's replicate `assemble_generic` logic.
    
    function make_tri_info(mesh, tri_idx)
        nodes = mesh.node
        elements = mesh.triangles
        tri_nodes = elements[:, tri_idx]
        v1 = nodes[:, tri_nodes[1]]
        v2 = nodes[:, tri_nodes[2]]
        v3 = nodes[:, tri_nodes[3]]
        vertices = hcat(v1, v2, v3)
        
        # Compute edges
        l1 = norm(v2 - v3)
        l2 = norm(v3 - v1)
        l3 = norm(v1 - v2)
        edgel = [l1, l2, l3]
        
        # Compute area and center
        center = (v1 + v2 + v3) / 3
        area = 0.5 * norm(cross(v2 - v1, v3 - v1))
        
        # Face normal
        n = cross(v2 - v1, v3 - v1)
        facen = n / norm(n)

        # Dummy values for others
        edgev = zeros(3, 3)
        edgen = zeros(3, 3)
        
        # Struct order: triID, tag, area, verticesID, vertices, center, facen̂, edgel, edgev̂, edgen̂, inBfsID, bfsSign
        return EMSuite.Geometry.TriangleInfo(
            tri_idx,                # triID
            0,                      # tag
            area,                   # area
            SVector{3, Int}(tri_nodes),  # verticesID
            SMatrix{3,3, Float64}(vertices), # vertices
            SVector{3, Float64}(center),     # center
            SVector{3, Float64}(facen),      # facen̂
            SVector{3, Float64}(edgel),      # edgel
            SMatrix{3,3, Float64}(edgev),    # edgev̂
            SMatrix{3,3, Float64}(edgen),    # edgen̂
            SVector{3, Int}([0,0,0]),    # inBfsID
            SVector{3, Int}([0,0,0])     # bfsSign
        )
    end
    
    tri1 = make_tri_info(mesh, tri_ids[1])
    tri2 = make_tri_info(mesh, tri_ids[2])
    
    # Calculate Components
    Z11_local = zeros(ComplexF64, 3, 3)
    Z22_local = zeros(ComplexF64, 3, 3)
    Z12_local = zeros(ComplexF64, 3, 3)
    Z21_local = zeros(ComplexF64, 3, 3)
    
    # Self Terms
    EMSuite.IntegralEquations.EFIEModule.calc_self_interaction!(Z11_local, efie, tri1)
    EMSuite.IntegralEquations.EFIEModule.calc_self_interaction!(Z22_local, efie, tri2)
    
    # Cross Terms
    EMSuite.IntegralEquations.EFIEModule.calc_near_interaction!(Z12_local, efie, tri1, tri2)
    EMSuite.IntegralEquations.EFIEModule.calc_near_interaction!(Z21_local, efie, tri2, tri1)
    
    # Apply Factor
    Z11_local .*= efie.factor
    Z22_local .*= efie.factor
    Z12_local .*= efie.factor
    Z21_local .*= efie.factor
    
    # Extract relevant local edges
    common_edge_len = bf_info["edge_length"]
    
    function get_local_edge_index(tri, len)
        for i in 1:3
            if abs(tri.edgel[i] - len) < 1e-6
                return i
            end
        end
        return 0
    end
    
    m1 = get_local_edge_index(tri1, common_edge_len)
    m2 = get_local_edge_index(tri2, common_edge_len)
    
    println("Local Edge Indices: m1=$m1, m2=$m2")
    
    z11 = Z11_local[m1, m1]
    z22 = Z22_local[m2, m2]
    z12 = Z12_local[m1, m2]
    z21 = Z21_local[m2, m1]
    
    println("Components (Raw):")
    println("  Z(T1, T1): ", z11)
    println("  Z(T2, T2): ", z22)
    println("  Z(T1, T2): ", z12)
    println("  Z(T2, T1): ", z21)
    
    z_total = z11 + z22 - z12 - z21
    println("Total Z[1,1] (Calculated): ", z_total)
    println("Legacy Z[1,1]: ", Z[1, 1])
    println("Difference: ", z_total - Z[1, 1])
    
end
