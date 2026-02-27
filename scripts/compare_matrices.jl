using JLD2
using LinearAlgebra
using Statistics
using Printf
using Plots

function compare_matrices()
    println("Loading matrices...")
    legacy = load(joinpath(@__DIR__, "..", "legacy_matrix.jld2"))
    emsuite = load(joinpath(@__DIR__, "..", "emsuite_matrix.jld2"))
    
    Z_leg = legacy["Z"]
    V_leg = legacy["V"]
    I_leg = legacy["I"]
    
    Z_em = emsuite["Z"]
    V_em = emsuite["V"]
    I_em = emsuite["I"]
    
    println("Legacy Z: $(size(Z_leg)), V: $(size(V_leg))")
    println("EMSuite Z: $(size(Z_em)), V: $(size(V_em))")
    
    # Compare Norms
    println("\n--- Norms ---")
    println("Legacy Z norm: $(norm(Z_leg))")
    println("EMSuite Z norm: $(norm(Z_em))")
    println("Ratio Z: $(norm(Z_em) / norm(Z_leg))")
    
    println("Legacy V norm: $(norm(V_leg))")
    println("EMSuite V norm: $(norm(V_em))")
    println("Ratio V: $(norm(V_em) / norm(V_leg))")
    
    println("Legacy I norm: $(norm(I_leg))")
    println("EMSuite I norm: $(norm(I_em))")
    println("Ratio I: $(norm(I_em) / norm(I_leg))")
    
    # Compare Elements
    println("\n--- Element Comparison ---")
    diff_Z = abs.(Z_em .- Z_leg)
    max_diff_Z = maximum(diff_Z)
    mean_diff_Z = mean(diff_Z)
    println("Max Diff Z: $max_diff_Z")
    println("Mean Diff Z: $mean_diff_Z")
    
    diff_V = abs.(V_em .- V_leg)
    max_diff_V = maximum(diff_V)
    mean_diff_V = mean(diff_V)
    println("Max Diff V: $max_diff_V")
    println("Mean Diff V: $mean_diff_V")
    
    # Check if it's just a scaling factor
    ratio_Z = Z_em ./ Z_leg
    mean_ratio_Z = mean(ratio_Z)
    std_ratio_Z = std(ratio_Z)
    println("Mean Ratio Z: $mean_ratio_Z")
    println("Std Ratio Z: $std_ratio_Z")
    
    ratio_V = V_em ./ V_leg
    mean_ratio_V = mean(ratio_V)
    std_ratio_V = std(ratio_V)
    println("Mean Ratio V: $mean_ratio_V")
    println("Std Ratio V: $std_ratio_V")
    
    # Check diagonal
    diag_leg = diag(Z_leg)
    diag_em = diag(Z_em)
    ratio_diag = diag_em ./ diag_leg
    println("Mean Ratio Diag: $(mean(ratio_diag))")
    
    println("\n--- Specific Elements ---")
    for i in 1:5
        for j in 1:5
            val_leg = Z_leg[i, j]
            val_em = Z_em[i, j]
            ratio = val_em / val_leg
            println("Z[$i,$j]: Leg=$(val_leg), EM=$(val_em), Ratio=$(ratio)")
        end
    end
    
    println("Legacy V[1]: $(V_leg[1])")
    println("EMSuite V[1]: $(V_em[1])")
    println("Ratio: $(V_em[1] / V_leg[1])")
    
    # Check if ordering is different
    # If ordering is different, the ratio will be random.
    # If ordering is same, ratio should be constant.
    
    if std_ratio_Z > 0.1
        println("WARNING: Z ratio is not constant. Ordering might be different.")
    end
    
end

compare_matrices()
