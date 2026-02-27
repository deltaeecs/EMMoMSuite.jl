using DelimitedFiles, Printf, Statistics

dir = joinpath("test_results", "emsuite_verification", "fullsphere")
b2 = readdlm(joinpath(dir, "B2_SMFIE_MLFMA.csv"), ',', Float64; skipstart=1)

# Try C3 first, then C1
for (label, fname) in [("C3_SCFIE_MLFMA", "C3_SCFIE_MLFMA.csv"), ("C1_SCFIE_Direct", "C1_SCFIE_Direct.csv")]
    fpath = joinpath(dir, fname)
    if isfile(fpath)
        ref = readdlm(fpath, ',', Float64; skipstart=1)
        println("B2 size: $(size(b2)), $label size: $(size(ref))")
        b2_e = b2[:, 3]
        ref_e = ref[:, 3]
        diff = b2_e .- ref_e
        rmse = sqrt(mean(diff .^ 2))
        @printf("B2 vs %s E-plane: Mean Diff=%.2f dB, RMSE=%.2f dB, Max|Diff|=%.2f dB\n", label, mean(diff), rmse, maximum(abs.(diff)))
        println("B2 range: $(round(minimum(b2_e),digits=1)) ~ $(round(maximum(b2_e),digits=1)) dBsm")
        println("$label range: $(round(minimum(ref_e),digits=1)) ~ $(round(maximum(ref_e),digits=1)) dBsm")
        println()
    end
end

println("Available files: ", readdir(dir))
