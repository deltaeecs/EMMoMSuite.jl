using Pkg

# Activate docs environment, make local package available, then resolve deps.
Pkg.activate(@__DIR__)
Pkg.develop(path = joinpath(@__DIR__, ".."))
Pkg.instantiate()
