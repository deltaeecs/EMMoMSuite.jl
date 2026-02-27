using Test
using EMSuite
using EMSuite.IO
using EMSuite.Geometry
using HDF5
using StaticArrays

@testset "IO" begin
    # Test save_RCS_txt
    theta = [0.0, pi/2]
    phi = [0.0, pi/2]
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
    
    # Test save_vtk
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
    
    # Test with data
    data_vec = [1.0]
    outfiles2 = save_vtk(filename_vtk, mesh, data_vec; data_name="TestVal")
    @test !isempty(outfiles2)
    for f in outfiles2
        rm(f)
    end
end
