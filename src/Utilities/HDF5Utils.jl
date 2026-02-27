module HDF5Utils

using HDF5
using SparseArrays

export save_sparse_matrix, load_sparse_matrix

function save_sparse_matrix(group::HDF5.Group, name::String, mat::SparseMatrixCSC)
    g = create_group(group, name)
    attributes(g)["type"] = "SparseMatrixCSC"
    g["m"] = mat.m
    g["n"] = mat.n
    g["colptr"] = mat.colptr
    g["rowval"] = mat.rowval
    g["nzval"] = mat.nzval
end

function load_sparse_matrix(group::HDF5.Group, name::String)
    g = group[name]
    if haskey(attributes(g), "type") && read(attributes(g)["type"]) == "SparseMatrixCSC"
        m = read(g["m"])
        n = read(g["n"])
        colptr = read(g["colptr"])
        rowval = read(g["rowval"])
        nzval = read(g["nzval"])
        return SparseMatrixCSC(m, n, colptr, rowval, nzval)
    else
        error("Group $name is not a SparseMatrixCSC")
    end
end

end
