module IO

include("ResultIO.jl")
include("VTKExport.jl")

using .ResultIO
using .VTKExport

export save_RCS_txt, save_results_hdf5, save_result
export save_vtk, save_vtk_multi

end
