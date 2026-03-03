module VTKExport

using WriteVTK
using StaticArrays
using ...CoreModule
using ...Geometry

export save_vtk, save_vtk_multi

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
function save_vtk(
    filename::String,
    mesh::TriangleMesh{IT,FT},
    data::Union{Nothing,AbstractVector{<:Real}} = nothing;
    data_name = "Data",
) where {IT,FT}

    # Prepare points
    points = mesh.node # 3 x Nv

    # Prepare cells
    # WriteVTK expects 1-based indexing for connectivity
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, mesh.triangles[:, i]) for i = 1:mesh.trinum]

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

function save_vtk(
    filename::String,
    mesh::TetrahedraMesh{IT,FT},
    data::Union{Nothing,AbstractVector{<:Real}} = nothing;
    data_name = "Data",
) where {IT,FT}

    points = mesh.node

    cells = [MeshCell(VTKCellTypes.VTK_TETRA, mesh.tetras[:, i]) for i = 1:mesh.tetnum]

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

# ─── HexahedraMesh ────────────────────────────────────────────────────────────

"""
    save_vtk(filename, mesh::HexahedraMesh, data=nothing; data_name="Data")

Export a hexahedral mesh (CHEXA convention, 8 nodes) to VTK unstructured grid
format. Node ordering follows the VTK_HEXAHEDRON convention (same as CHEXA).
"""
function save_vtk(
    filename::String,
    mesh::HexahedraMesh{IT,FT},
    data::Union{Nothing,AbstractVector{<:Real}} = nothing;
    data_name = "Data",
) where {IT,FT}
    points = mesh.node
    cells  = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, mesh.hexes[:, i]) for i = 1:mesh.hexnum]
    vtk    = vtk_grid(filename, points, cells)

    if data !== nothing
        if length(data) == size(points, 2)
            vtk_point_data(vtk, data, data_name)
        elseif length(data) == mesh.hexnum
            vtk_cell_data(vtk, data, data_name)
        else
            @warn "Data length mismatch. Skipping data."
        end
    end

    outfiles = vtk_save(vtk)
    return outfiles
end

# ─── Complex vector field ─────────────────────────────────────────────────────

"""
    save_vtk(filename, mesh, vector_data; field_name="E_field", save_mode=:real_imag)

Export a complex 3-component vector field dataset to VTK.

`vector_data` must be a `Vector{SVector{3,Complex{FT}}}` aligned with **nodes**
(point data) or **elements** (cell data) of `mesh`.

`save_mode`:
- `:real_imag` (default) — saves `<field_name>_Re` and `<field_name>_Im`
- `:magnitude` — saves `<field_name>_Abs` and `<field_name>_Phase_deg`
"""
function save_vtk(
    filename::String,
    mesh::AbstractMesh,
    vector_data::AbstractVector{SVector{3,CT}};
    field_name::String = "E_field",
    save_mode::Symbol  = :real_imag,
) where {CT<:Complex}
    # Determine whether it's triangles/tets/hexes and get the basic vtk grid
    # We delegate by calling the fundamental save_vtk(no data), then add fields.
    # Since WriteVTK doesn't allow re-opening, we build from scratch:
    nv = size(mesh.node, 2)

    if mesh isa TriangleMesh
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, mesh.triangles[:, i]) for i = 1:mesh.trinum]
    elseif mesh isa TetrahedraMesh
        cells = [MeshCell(VTKCellTypes.VTK_TETRA, mesh.tetras[:, i]) for i = 1:mesh.tetnum]
    elseif mesh isa HexahedraMesh
        cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, mesh.hexes[:, i]) for i = 1:mesh.hexnum]
    else
        error("save_vtk: unsupported mesh type $(typeof(mesh)) for complex vector export")
    end

    vtk = vtk_grid(filename, mesh.node, cells)

    # Attach vector data as 3×N real arrays
    ET = eltype(eltype(vector_data))
    FT = ET <: Complex ? real(ET) : ET
    n  = length(vector_data)

    if save_mode == :real_imag
        data_re  = Matrix{FT}(undef, 3, n)
        data_im  = Matrix{FT}(undef, 3, n)
        for i in 1:n
            data_re[:,i] .= real.(vector_data[i])
            data_im[:,i] .= imag.(vector_data[i])
        end
        if n == nv
            vtk_point_data(vtk, data_re, field_name * "_Re")
            vtk_point_data(vtk, data_im, field_name * "_Im")
        else
            vtk_cell_data(vtk, data_re, field_name * "_Re")
            vtk_cell_data(vtk, data_im, field_name * "_Im")
        end
    elseif save_mode == :magnitude
        data_abs   = Matrix{FT}(undef, 3, n)
        data_phase = Matrix{FT}(undef, 3, n)
        for i in 1:n
            for k in 1:3
                data_abs[k,i]   = abs(vector_data[i][k])
                data_phase[k,i] = rad2deg(angle(vector_data[i][k]))
            end
        end
        if n == nv
            vtk_point_data(vtk, data_abs,   field_name * "_Abs")
            vtk_point_data(vtk, data_phase, field_name * "_Phase_deg")
        else
            vtk_cell_data(vtk, data_abs,   field_name * "_Abs")
            vtk_cell_data(vtk, data_phase, field_name * "_Phase_deg")
        end
    else
        error("save_vtk: unknown save_mode=$save_mode (use :real_imag or :magnitude)")
    end

    return vtk_save(vtk)
end

# ─── Multi-field export ───────────────────────────────────────────────────────

"""
    save_vtk_multi(filename, mesh; point_data=Dict(), cell_data=Dict())

Export a mesh together with multiple scalar or vector data sets in a single
`.vtu` file.  Keys become the field names in ParaView.

Each value in `point_data` / `cell_data` must be a scalar `Vector{<:Real}` or
a 3D vector `Vector{SVector{3}}` (or a `3×N Matrix`).
"""
function save_vtk_multi(
    filename::String,
    mesh::AbstractMesh;
    point_data::Dict{String} = Dict{String,Any}(),
    cell_data::Dict{String}  = Dict{String,Any}(),
)
    nv = size(mesh.node, 2)

    if mesh isa TriangleMesh
        cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, mesh.triangles[:, i]) for i = 1:mesh.trinum]
    elseif mesh isa TetrahedraMesh
        cells = [MeshCell(VTKCellTypes.VTK_TETRA, mesh.tetras[:, i]) for i = 1:mesh.tetnum]
    elseif mesh isa HexahedraMesh
        cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, mesh.hexes[:, i]) for i = 1:mesh.hexnum]
    else
        error("save_vtk_multi: unsupported mesh type $(typeof(mesh))")
    end

    vtk = vtk_grid(filename, mesh.node, cells)

    for (name, val) in point_data
        if val isa AbstractVector{<:SVector}
            mat = reduce(hcat, val)   # 3×N or similar
            vtk_point_data(vtk, mat, name)
        else
            vtk_point_data(vtk, vec(val), name)
        end
    end

    for (name, val) in cell_data
        if val isa AbstractVector{<:SVector}
            mat = reduce(hcat, val)
            vtk_cell_data(vtk, mat, name)
        else
            vtk_cell_data(vtk, vec(val), name)
        end
    end

    return vtk_save(vtk)
end

end
