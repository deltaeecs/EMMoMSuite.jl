using EMSuite
using LinearAlgebra
using StaticArrays
using Test
using Printf
using SparseArrays

function run_benchmark()
    println("==================================================")
    println("   Benchmark: MLFMA vs Direct Solver (Accuracy)   ")
    println("==================================================")

    # 1. Parameters
    freq = 300e6 # 300 MHz
    lambda = 299792458.0 / freq
    radius = 0.5 # Small radius to keep direct solver fast
    
    println(@sprintf("Frequency: %.2f MHz", freq/1e6))
    println(@sprintf("Wavelength: %.4f m", lambda))
    println(@sprintf("Radius: %.4f m", radius))

    # 2. Mesh Generation
    # Use a coarse mesh
    n_theta = 12
    n_phi = 24
    
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    println(@sprintf("Vertices: %d", num_vertices(mesh)))
    println(@sprintf("Elements: %d", num_elements(mesh)))

    # 3. Basis Setup
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println(@sprintf("Unknowns: %d", N))

    # 4. Direct Solver (EFIE)
    println("\n[1/2] Running Direct Solver...")
    efie = EFIE(freq)
    
    t_direct = @elapsed begin
        Z_direct = assemble_impedance_matrix(efie, basis)
    end
    println(@sprintf("      Assembly: %.4f s", t_direct))
    
    source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
    V = excitation_vector(efie, source, basis)
    
    # Solve Direct
    I_direct = Z_direct \ V
    println("      Solved Direct.")

    # 5. MLFMA Solver
    println("\n[2/2] Running MLFMA Solver...")
    
    # Leaf box size: 0.25 lambda
    leaf_size = 0.25 * lambda
    
    t_mlfma_setup = @elapsed begin
        Z_mlfma = MLFMAOperator(efie, basis, leaf_size)
    end
    println(@sprintf("      MLFMA Setup: %.4f s", t_mlfma_setup))
    
    # --- DEBUG: Check Near Field Consistency ---
    println("\n[DEBUG] Checking Near Field Consistency...")
    sorted_ids = Z_mlfma.sorted_ids
    println("Debug: sorted_ids[1] = ", sorted_ids[1])
    # Z_near is in Original Basis (based on MLFMAOperator implementation)
    # So we compare directly with Z_direct
    Z_near = Z_mlfma.Z_near
    
    # Find indices where Z_near is non-zero
    rows, cols, vals = findnz(Z_near)
    
    max_diff = 0.0
    sum_diff = 0.0
    count = 0
    
    for i in 1:length(vals)
        r = rows[i]
        c = cols[i]
        val_near = vals[i]
        val_direct = Z_direct[r, c]
        
        diff = abs(val_near - val_direct)
        if diff > max_diff
            max_diff = diff
        end
        
        if diff > 1.0 && count < 5
             println(@sprintf("      Mismatch at (%d, %d): Near=%.2e+j%.2e, Direct=%.2e+j%.2e, Diff=%.2e", 
                r, c, real(val_near), imag(val_near), real(val_direct), imag(val_direct), diff))
        end
        
        sum_diff += diff
        count += 1
    end
    
    println(@sprintf("      Checked %d near-field entries.", count))
    println(@sprintf("      Max Difference: %.2e", max_diff))
    println(@sprintf("      Avg Difference: %.2e", sum_diff / count))
    
    if max_diff > 1e-10
        println("      WARNING: Near field mismatch detected!")
    else
        println("      SUCCESS: Near field matches Direct Solver.")
    end
    # -------------------------------------------

    # Solve MLFMA (GMRES)
    # Use ILU Preconditioner to ensure convergence
    println("      Computing ILU Preconditioner...")
    P = ILUPreconditioner(Z_mlfma.Z_near, τ=0.01)
    
    solver = GMRESSolver(tol=1e-5, maxiter=200, verbose=false)
    
    t_mlfma_solve = @elapsed begin
        I_mlfma = solve!(solver, Z_mlfma, V, Pl=P)
    end
    println(@sprintf("      Solved MLFMA in %.4f s", t_mlfma_solve))
    
    # 6. Comparison
    norm_direct = norm(I_direct)
    norm_mlfma = norm(I_mlfma)
    println(@sprintf("      Norm Direct: %.4e", norm_direct))
    println(@sprintf("      Norm MLFMA:  %.4e", norm_mlfma))
    
    # Check Residual
    res = V - Z_mlfma * I_mlfma
    norm_res = norm(res) / norm(V)
    println(@sprintf("      Relative Residual: %.4e", norm_res))
    
    diff = norm(I_mlfma - I_direct) / norm(I_direct)
    println("\n--------------------------------------------------")
    println(@sprintf("Relative Error (Norm): %.2e", diff))
    println("--------------------------------------------------")
    
    if diff < 0.1
        println("SUCCESS: MLFMA matches Direct Solver (Error < 10%).")
    else
        println("WARNING: High discrepancy between MLFMA and Direct Solver.")
    end
end

run_benchmark()
