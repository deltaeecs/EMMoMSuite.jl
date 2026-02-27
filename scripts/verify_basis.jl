using EMSuite
using JLD2
using StaticArrays
using LinearAlgebra
using Printf

# Load reference data
# Note: We need to handle the types from MoM_Basics which are not available here unless we define them or use a generic loader.
# JLD2 might complain if types are missing.
# We can try to load into a Dict or use @load.
# If MoM_Basics types are not available, JLD2 will load them as ReconstructedTypes.

# We need to define a minimal structure to match if JLD2 fails, but let's try loading first.
# To avoid type issues, we can just load the raw data (nodes, triangles) and the basis info as "Any" or reconstructed.

data = load("../reference_data.jld2")
nodes = data["nodes"]
triangles = data["triangles"]
ref_stbfsInfo = data["stbfsInfo"]
ref_trisInfo = data["trisInfo"]

println("Loaded reference data.")
println("Reference RWG count: ", length(ref_stbfsInfo))

# Create EMSuite Mesh
mesh = TriangleMesh(size(triangles, 2), nodes, triangles)

# Create EMSuite Basis
basis = RWGBasis(mesh)
println("EMSuite RWG count: ", num_basis(basis))

if length(ref_stbfsInfo) != num_basis(basis)
    println("ERROR: Number of basis functions mismatch!")
    exit(1)
end

# Compare
# Build map for reference basis: Center -> Index
ref_map = Dict{SVector{3, Float64}, Int}()

# Helper to extract center from ReconstructedType or Struct
function get_center(bf)
    # If it's a struct
    if hasproperty(bf, :center)
        return bf.center
    end
    # If it's a JLD2 ReconstructedType, we might need to access fields differently
    # But usually property access works if fields match.
    # Let's assume property access works.
    return bf.center
end

for (i, bf) in enumerate(ref_stbfsInfo)
    c = get_center(bf)
    # Round to avoid float issues
    c_rounded = round.(c, digits=6)
    ref_map[c_rounded] = i
end

matched_count = 0
mismatch_count = 0

for i in 1:num_basis(basis)
    bf = basis.functions[i]
    c = bf.center
    c_rounded = round.(c, digits=6)
    
    if haskey(ref_map, c_rounded)
        ref_idx = ref_map[c_rounded]
        ref_bf = ref_stbfsInfo[ref_idx]
        
        if ref_idx != i
            println("Index mismatch for BF $i: Ref index is $ref_idx")
            global mismatch_count += 1
        end
        
        # Compare Edge Length
        # Access edgel from ref_bf
        ref_len = ref_bf.edgel
        if abs(bf.edge_length - ref_len) > 1e-6
            println("Mismatch in Edge Length for BF $i (Ref $ref_idx): EMSuite=$(bf.edge_length), Ref=$(ref_len)")
            global mismatch_count += 1
        end
        
        # Compare Support
        # Ref: inGeo (Vector of 2 Ints)
        ref_support = ref_bf.inGeo
        
        # EMSuite: support (SVector of 2 Ints)
        # Check for orientation flip
        is_flipped = false
        if bf.support[1] == ref_support[1] && bf.support[2] == ref_support[2]
            # Match
        elseif bf.support[1] == ref_support[2] && bf.support[2] == ref_support[1]
            # Flipped
            is_flipped = true
            println("Orientation Mismatch for BF $i (Ref $ref_idx): EMSuite=$(bf.support), Ref=$(ref_support)")
            global mismatch_count += 1
        else
            println("Mismatch in Support for BF $i (Ref $ref_idx): EMSuite=$(bf.support), Ref=$(ref_support)")
            global mismatch_count += 1
        end
        
        global matched_count += 1
    else
        println("Could not find match for BF $i at center $c")
        global mismatch_count += 1
    end
end

println("Verification Complete.")
println("Matched: $matched_count")
println("Mismatches: $mismatch_count")

if mismatch_count == 0
    println("BASIS FUNCTIONS VERIFIED: SUCCESS")
else
    println("BASIS FUNCTIONS VERIFIED: FAILURE")
end
