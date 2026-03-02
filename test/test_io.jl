using Test
using EMSuite
using EMSuite.IO
using EMSuite.Geometry
using EMSuite.Utilities: save_sparse_matrix, load_sparse_matrix
using HDF5
using SparseArrays
using StaticArrays

@testset "IO" begin
    # Test save_RCS_txt
    theta = [0.0, pi / 2]
    phi = [0.0, pi / 2]
    rcs_data = [1.0 2.0; 3.0 4.0]

    filename_txt = "test_rcs.txt"
    save_RCS_txt(filename_txt, theta, phi, rcs_data)

    @test isfile(filename_txt)
    rm(filename_txt)

    # Test save_results_hdf5
    filename_h5 = "test_results.h5"
    data = Dict("a" => 1, "b" => [1, 2, 3])
    save_results_hdf5(filename_h5, data)

    @test isfile(filename_h5)

    # Load back to verify
    h5open(filename_h5, "r") do file
        @test read(file["a"]) == 1
        @test read(file["b"]) == [1, 2, 3]
    end

    rm(filename_h5)

    # Test save_results_hdf5 kwargs version
    filename_h5b = "test_results_kw.h5"
    save_results_hdf5(filename_h5b; x=42, y=[7.0, 8.0])
    @test isfile(filename_h5b)
    h5open(filename_h5b, "r") do file
        @test read(file["x"]) == 42
        @test read(file["y"]) ≈ [7.0, 8.0]
    end
    rm(filename_h5b)

    # Test save_vtk (TriangleMesh)
    v1 = SVector(0.0, 0.0, 0.0)
    v2 = SVector(1.0, 0.0, 0.0)
    v3 = SVector(0.0, 1.0, 0.0)
    node = Matrix(hcat(v1, v2, v3))
    triangles = reshape([1, 2, 3], 3, 1)
    mesh = TriangleMesh(1, node, triangles)

    filename_vtk = "test_mesh" # WriteVTK adds .vtu
    outfiles = save_vtk(filename_vtk, mesh)

    @test !isempty(outfiles)
    @test isfile(outfiles[1])

    # Clean up
    for f in outfiles
        rm(f)
    end

    # Test with data (point data, same length as vertices)
    data_vec = [1.0, 2.0, 3.0]  # 3 point values
    outfiles2 = save_vtk(filename_vtk, mesh, data_vec; data_name = "PointVal")
    @test !isempty(outfiles2)
    for f in outfiles2
        rm(f)
    end

    # Test with cell data (same length as elements)
    data_cell = [99.0]  # 1 cell value
    outfiles3 = save_vtk(filename_vtk, mesh, data_cell; data_name = "CellVal")
    @test !isempty(outfiles3)
    for f in outfiles3
        rm(f)
    end

    # Test mismatched data (should warn and skip)
    data_bad = [1.0, 2.0]  # neither 3 vertices nor 1 cell
    outfiles4 = save_vtk(filename_vtk, mesh, data_bad; data_name = "BadData")
    @test !isempty(outfiles4)  # Still saves, just without data
    for f in outfiles4
        rm(f)
    end

    # Test save_vtk (TetrahedraMesh)
    nodes_tet = [
        0.0  1.0  0.0  0.0;
        0.0  0.0  1.0  0.0;
        0.0  0.0  0.0  1.0
    ]
    tetras = reshape([1, 2, 3, 4], 4, 1)
    mesh_tet = TetrahedraMesh(1, nodes_tet, tetras)

    outfiles_tet = save_vtk(filename_vtk, mesh_tet)
    @test !isempty(outfiles_tet)
    for f in outfiles_tet
        rm(f)
    end

    # TetrahedraMesh with cell data
    outfiles_tet2 = save_vtk(filename_vtk, mesh_tet, [42.0]; data_name = "TetCell")
    @test !isempty(outfiles_tet2)
    for f in outfiles_tet2
        rm(f)
    end

    # TetrahedraMesh with mismatched data (should warn)
    outfiles_tet3 = save_vtk(filename_vtk, mesh_tet, [1.0, 2.0]; data_name = "TetBad")
    @test !isempty(outfiles_tet3)
    for f in outfiles_tet3
        rm(f)
    end

    # Test HDF5Utils: save_sparse_matrix / load_sparse_matrix
    filename_sparse = "test_sparse.h5"
    # Real sparse matrix — must pass HDF5.Group, not File; create a container group
    A_real = sparse([1, 2, 3, 1], [1, 2, 3, 3], [1.0, 2.5, 3.0, 4.0], 4, 4)
    h5open(filename_sparse, "w") do file
        grp = create_group(file, "matrices")
        save_sparse_matrix(grp, "A_real", A_real)
    end
    @test isfile(filename_sparse)
    A_loaded = h5open(filename_sparse, "r") do file
        load_sparse_matrix(file["matrices"], "A_real")
    end
    @test A_loaded ≈ A_real
    @test size(A_loaded) == size(A_real)
    @test nnz(A_loaded) == nnz(A_real)
    rm(filename_sparse)

    # Complex sparse matrix
    filename_sparse2 = "test_sparse_complex.h5"
    B_complex = sparse([1, 2], [1, 2], [1.0 + 2.0im, 3.0 - 1.0im], 3, 3)
    h5open(filename_sparse2, "w") do file
        grp = create_group(file, "mats")
        save_sparse_matrix(grp, "B_complex", B_complex)
    end
    @test isfile(filename_sparse2)
    B_loaded = h5open(filename_sparse2, "r") do file
        load_sparse_matrix(file["mats"], "B_complex")
    end
    @test B_loaded ≈ B_complex
    rm(filename_sparse2)

    # load_sparse_matrix: type mismatch error (no "type" attribute)
    filename_bad = "test_sparse_bad.h5"
    h5open(filename_bad, "w") do file
        grp = create_group(file, "bad")
        g = create_group(grp, "NotAMatrix")
        g["values"] = [1, 2, 3]
        # No "type" attribute set → should throw ErrorException
    end
    @test_throws ErrorException h5open(filename_bad, "r") do file
        load_sparse_matrix(file["bad"], "NotAMatrix")
    end
    rm(filename_bad)

end
