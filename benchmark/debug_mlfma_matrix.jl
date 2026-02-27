using EMSuite
using LinearAlgebra
using SparseArrays
using Printf
using Statistics

# Paths
const MOM_ALLINONE_DIR = "f:/OneDrive/MoM/MoM_AllinOne"

function debug_mlfma_matrix()
    println("=== Debugging MLFMA Matrix Elements ===")
    
    # 1. Load Mesh (Small mesh for debugging)
    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles/jet_100MHz.nas")
    println("Loading mesh: $mesh_file")
    
    mesh = read_nas_mesh(mesh_file, scale=1.0)
    
    # 2. Parameters
    freq = 1e8
    set_frequency!(freq)
    lambda = 299792458.0 / freq
    k = 2π / lambda
    
    # 3. Basis
    println("Setting up RWG basis...")
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Number of unknowns: $N")
    
    # 4. Equation
    println("Setting up EFIE...")
    efie = EFIE(freq)
    
    # 5. MLFMA Operator
    println("Setting up MLFMA Operator...")
    leaf_size = 0.25 * lambda
    Z_mlfma_op = MLFMAOperator(efie, basis, leaf_size)
    
    # 6. Direct Matrix (Reference)
    println("Assembling Direct Matrix (Reference)...")
    t_start = time()
    Z_direct = assemble_impedance_matrix(efie, basis)
    println("Direct assembly took $(time() - t_start) seconds.")
    
    # 7. Compute MLFMA Matrix Column by Column (Sample)
    println("Computing MLFMA Matrix via matvecs (Sampling 50 columns)...")
    
    # Sample random columns
    n_samples = 50
    sample_indices = sort(rand(1:N, n_samples))
    
    Z_mlfma_dense = zeros(ComplexF64, N, n_samples)
    Z_direct_sample = Z_direct[:, sample_indices]
    
    # Use threads for speed
    Threads.@threads for k in 1:n_samples
        j = sample_indices[k]
        x = zeros(ComplexF64, N)
        x[j] = 1.0
        y = Z_mlfma_op * x
        Z_mlfma_dense[:, k] = y
        if k % 10 == 0
            println("Computed column $k / $n_samples")
        end
    end
    
    # 8. Compare
    println("\n=== Comparison Results ===")
    
    # Element-wise difference
    Diff = Z_mlfma_dense .- Z_direct_sample
    
    # Separate into Near and Far
    rows, cols, vals = findnz(Z_mlfma_op.Z_near)
    near_indices = Set(zip(rows, cols))
    
    diff_near = ComplexF64[]
    diff_far = ComplexF64[]
    
    val_direct_far = ComplexF64[]
    val_mlfma_far = ComplexF64[]
    
    for k in 1:n_samples
        j = sample_indices[k]
        for i in 1:N
            if (i, j) in near_indices
                push!(diff_near, Diff[i, k])
            else
                push!(diff_far, Diff[i, k])
                push!(val_direct_far, Z_direct_sample[i, k])
                push!(val_mlfma_far, Z_mlfma_dense[i, k])
            end
        end
    end
    
    println("Near Field Terms: $(length(diff_near))")
    println("Far Field Terms: $(length(diff_far))")
    
    if !isempty(diff_near)
        max_diff_near = maximum(abs.(diff_near))
        rmse_near = sqrt(mean(abs2.(diff_near)))
        println("Max Diff Near: $max_diff_near")
        println("RMSE Near: $rmse_near")
    end
    
    if !isempty(diff_far)
        max_diff_far = maximum(abs.(diff_far))
        rmse_far = sqrt(mean(abs2.(diff_far)))
        println("Max Diff Far: $max_diff_far")
        println("RMSE Far: $rmse_far")
        
        # Check Ratio for Far Field
        mask = abs.(val_direct_far) .> 1e-10
        if any(mask)
            ratios = abs.(val_mlfma_far[mask]) ./ abs.(val_direct_far[mask])
            mean_ratio = mean(ratios)
            println("Mean Ratio (MLFMA / Direct) for Far Field: $mean_ratio")
            println("Expected Ratio if 4pi missing: $(1/(4π)) = $(1/(4π))")
            println("Expected Ratio if 4pi extra: $(4π) = $(4π)")
        end
    end
    
    # Save a sample of comparison
    open("matrix_debug_sample.txt", "w") do io
        println(io, "Row, Col, Direct_Abs, MLFMA_Abs, Ratio, Type")
        for _ in 1:50
            k = rand(1:n_samples)
            j = sample_indices[k]
            i = rand(1:N)
            type = (i, j) in near_indices ? "NEAR" : "FAR"
            v_dir = Z_direct_sample[i, k]
            v_mlfma = Z_mlfma_dense[i, k]
            ratio = abs(v_dir) > 1e-12 ? abs(v_mlfma) / abs(v_dir) : 0.0
            println(io, "$i, $j, $(abs(v_dir)), $(abs(v_mlfma)), $ratio, $type")
        end
    end
    println("Saved sample to matrix_debug_sample.txt")
end

debug_mlfma_matrix()
