using EMSuite.CoreModule
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.FastAlgorithms.MLFMA
using LinearAlgebra
using StaticArrays
using SparseArrays

function verify_matrix_elements()
    println("=== Verifying Matrix Elements ===")
    
    # 1. Setup
    # Use a slightly larger mesh to ensure we have both near and far interactions
    # Sphere radius 1.0
    # Frequency 300 MHz (lambda = 1.0)
    # Box size 0.5 -> 2 levels
    
    freq = 300e6
    k = 2π * freq / 299792458.0
    operator = EFIE(freq)
    
    # Create a simple mesh (Octahedron)
    nodes = [
        1.0 0.0 0.0;
        -1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 -1.0 0.0;
        0.0 0.0 1.0;
        0.0 0.0 -1.0
    ]
    
    # 8 faces
    faces = [
        1 3 5; 1 5 4; 1 4 6; 1 6 3;
        2 5 3; 2 4 5; 2 6 4; 2 3 6
    ]
    
    nt = size(faces, 1)
    mesh = TriangleMesh(nt, collect(nodes'), collect(faces'))
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("Number of unknowns: ", N)
    
    # 2. Dense Matrix
    println("Computing Dense Matrix...")
    Z_dense = assemble_impedance_matrix(operator, basis)
    
    # 3. MLFMA Operator
    println("Setting up MLFMA...")
    # Force small box size to ensure multiple levels
    # Extent is 2.0. Box size 0.4 should give ~5 boxes -> Level 3?
    mlfma = MLFMAOperator(operator, basis, 0.4)
    
    # 4. Compare Near Field
    println("\n--- Comparing Near Field ---")
    Z_near = mlfma.Z_near
    # Z_near is in the ORIGINAL basis, so we don't need sorted_ids to map indices
    
    max_diff_near = 0.0
    count_near = 0
    
    rows, cols, vals = findnz(Z_near)
    for k in 1:length(vals)
        r_orig = rows[k]
        c_orig = cols[k]
        val_near = vals[k]
        
        val_dense = Z_dense[r_orig, c_orig]
        
        diff = abs(val_near - val_dense)
        if diff > max_diff_near
            max_diff_near = diff
        end
        count_near += 1
        
        if diff > 1e-10
            println("Mismatch at ($r_orig, $c_orig): Near=$val_near, Dense=$val_dense, Diff=$diff")
        end
    end
    
    println("Checked $count_near near-field elements.")
    println("Max Near Field Difference: $max_diff_near")
    
    # 5. Compare Far Field (Column by Column)
    println("\n--- Comparing Far Field ---")
    # Pick a source basis function
    src_idx_orig = 1
    println("Testing source basis function: $src_idx_orig")
    
    x = zeros(ComplexF64, N)
    x[src_idx_orig] = 1.0
    
    # Dense result
    y_dense = Z_dense * x
    
    # MLFMA result
    y_mlfma = zeros(ComplexF64, N)
    mul!(y_mlfma, mlfma, x)
    
    # Separate Near and Far in Dense
    # We need to know which entries are "near" for this source
    # We can use the Z_near structure
    
    col_near = Z_near[:, src_idx_orig] # Sparse vector
    
    y_near_dense = zeros(ComplexF64, N)
    
    # Map back to original indices
    nz_rows, nz_vals = findnz(col_near)
    for k in 1:length(nz_rows)
        r_orig = nz_rows[k]
        val = nz_vals[k]
        y_near_dense[r_orig] = val
    end
    
    y_far_dense = y_dense - y_near_dense
    
    # MLFMA Far part
    # We can compute it by subtracting the near part from total MLFMA
    # Or we can instrument MLFMA to return only far part.
    # For now, let's assume y_mlfma = y_near_mlfma + y_far_mlfma
    # Since we verified Near Field matches (hopefully), 
    # y_far_mlfma ≈ y_mlfma - y_near_dense
    
    y_far_mlfma = y_mlfma - y_near_dense
    
    println("Norm Dense Total: ", norm(y_dense))
    println("Norm MLFMA Total: ", norm(y_mlfma))
    println("Total Difference: ", norm(y_dense - y_mlfma))
    println("Norm Far Dense:   ", norm(y_far_dense))
    println("Norm Far MLFMA:   ", norm(y_far_mlfma))
    
    diff_far = norm(y_far_dense - y_far_mlfma)
    println("Far Field Difference Norm: ", diff_far)
    
    if norm(y_far_dense) > 1e-12
        println("Relative Far Error: ", diff_far / norm(y_far_dense))
        println("Ratio Far (MLFMA/Dense): ", norm(y_far_mlfma) / norm(y_far_dense))
    end
    
    # Detailed Far Field Comparison
    println("\nDetailed Far Field Elements:")
    for i in 1:N
        val_d = y_far_dense[i]
        val_m = y_far_mlfma[i]
        if abs(val_d) > 1e-12 || abs(val_m) > 1e-12
            println("Obs $i: Dense=$val_d, MLFMA=$val_m, Ratio=$(abs(val_m)/abs(val_d)), PhaseDiff=$(angle(val_m/val_d))")
        end
    end

end

verify_matrix_elements()
