# test_cov_io.jl — IO 模块边界覆盖（VTK/ResultIO 变体）
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.IO
using EMMoMSuite.CoreModule
using LinearAlgebra

@testset "VTK 导出变体" begin
    mesh = generate_sphere_mesh(0.5, 4, 8)
    path = joinpath(tempdir(), "cov_mesh_$(getpid()).vtu")
    save_vtk(path, mesh)
    @test isfile(path)
    rm(path; force = true)
    # 四面体网格 VTK
    vol = generate_box_tet_mesh(0.5, 0.5, 0.5, 1, 1, 1)
    path2 = joinpath(tempdir(), "cov_tet_$(getpid()).vtu")
    save_vtk(path2, vol)
    @test isfile(path2)
    rm(path2; force = true)
end

@testset "ResultIO" begin
    # RCS CSV/TXT
    path = joinpath(tempdir(), "cov_rcs_$(getpid())")
    freqs = [1e9, 2e9]
    theta = [0.0, pi / 2]
    phi = [0.0]
    comp = rand(2, length(theta), length(phi))
    total = rand(length(theta), length(phi))
    dB = 10 * log10.(total .+ 1e-30)
    save_RCS_csv(path, theta, phi, comp, total, dB)
    @test isfile(path * ".csv")
    rm(path * ".csv"; force = true)
    path2 = joinpath(tempdir(), "cov_rcs_$(getpid()).txt")
    save_RCS_txt(path2, theta, phi, rand(length(theta), length(phi)))
    @test isfile(path2)
    rm(path2; force = true)
    # HDF5 结果保存
    path3 = joinpath(tempdir(), "cov_res_$(getpid()).h5")
    save_results_hdf5(path3, Dict("I" => rand(ComplexF64, 3), "freq" => 1e9))
    @test isfile(path3)
    rm(path3; force = true)
end
