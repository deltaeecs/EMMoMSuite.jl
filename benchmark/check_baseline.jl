# Quick look at Legacy baseline data
using CSV, DataFrames, Statistics

df = CSV.read(joinpath(@__DIR__, "..", "test_results", "legacy_baseline", "SEFIE_Direct_Jet.csv"), DataFrame)
println("Columns: ", names(df))
println("Size: ", size(df))
println("\nFirst 5 rows:")
println(first(df, 5))
println("\nTheta range: ", df.Theta_Rad[1], " to ", df.Theta_Rad[end])
println("Phi=0 dBsm: min=$(minimum(df.RCS_Phi0_dBsm)), max=$(maximum(df.RCS_Phi0_dBsm))")
println("Phi=90 dBsm: min=$(minimum(df.RCS_Phi90_dBsm)), max=$(maximum(df.RCS_Phi90_dBsm))")

# Key angles
for i in 1:length(df.Theta_Rad)
    td = round(df.Theta_Rad[i] * 180 / π, digits=1)
    if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
        println("θ=$(td)°: Phi=0 → $(df.RCS_Phi0_dBsm[i]) dBsm, Phi=90 → $(df.RCS_Phi90_dBsm[i]) dBsm")
    end
end
