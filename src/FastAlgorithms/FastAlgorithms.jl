module FastAlgorithms

include("MLFMA/MLFMA.jl")
include("Lebedev/Lebedev.jl")

using .MLFMA
using .Lebedev

export MLFMA, Lebedev, MLFMAOperator

end
