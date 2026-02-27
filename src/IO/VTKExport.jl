module VTKExport

using WriteVTK
using StaticArrays
using ...Geometry

export save_vtk

"""
    save_vtk(filename, mesh, data=nothing; data_name="Data")

Export mesh and associated data to VTK format (.vtu) for visualization in ParaView or VisIt.

# Arguments
- `filename`: Output filename (without extension).
- `mesh`: The mesh object (`TriangleMesh` or `TetrahedraMesh`).
- `data`: Optional data vector to attach to the mesh.
- `data_name`: Name of the data field in the VTK file.

# Data Handling
- If `length(data) == num_vertices(mesh)`, it is treated as **Point Data** (interpolated across elements).
- If `length(data) == num_elements(mesh)`, it is treated as **Cell Data** (constant per element).
- If dimensions mismatch, a warning is issued and data is skipped.

# Returns
- `outfiles`: List of generated files.
"""
function save_vtk(filename::String, mesh::TriangleMesh{IT, FT}, data::Union{Nothing, AbstractVector}=nothing; data_name="Data") where {IT, FT}
    
    # Prepare points
    points = mesh.node # 3 x Nv
    
    # Prepare cells
    # WriteVTK expects 1-based indexing for connectivity
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, mesh.triangles[:, i]) for i in 1:mesh.trinum]
    
    # Create grid
    vtk = vtk_grid(filename, points, cells)
    
    # Add data if present
    if data !== nothing
        # Check if data is point data or cell data
        if length(data) == size(points, 2)
            vtk_point_data(vtk, data, data_name)
        elseif length(data) == mesh.trinum
            vtk_cell_data(vtk, data, data_name)
        else
            @warn "Data length ($(length(data))) does not match number of points ($(size(points, 2))) or cells ($(mesh.trinum)). Skipping data."
        end
    end
    
    # Save
    outfiles = vtk_save(vtk)
    return outfiles
end

function save_vtk(filename::String, mesh::TetrahedraMesh{IT, FT}, data::Union{Nothing, AbstractVector}=nothing; data_name="Data") where {IT, FT}
    
    points = mesh.node
    
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, mesh.tetras[:, i]) for i in 1:mesh.tetnum]
    
    vtk = vtk_grid(filename, points, cells)
    
    if data !== nothing
        if length(data) == size(points, 2)
            vtk_point_data(vtk, data, data_name)
        elseif length(data) == mesh.tetnum
            vtk_cell_data(vtk, data, data_name)
        else
            @warn "Data length mismatch. Skipping data."
        end
    end
    
    outfiles = vtk_save(vtk)
    return outfiles
end

end
