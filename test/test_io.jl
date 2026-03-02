using Test
using EMSuite
using EMSuite.IO
using EMSuite.Geometry
using HDF5
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
end
