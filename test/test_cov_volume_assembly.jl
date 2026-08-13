# test_cov_volume_assembly.jl — VolumeAssembly.jl 轻量 P=1 MPI 装配覆盖
#
# 覆盖 src/Parallel/MPI/VolumeAssembly.jl 的全部 public 装配方法：
#   VEFIE + SWGBasis/PWCBasis/PWCHexBasis
#   SCFIE + RWGBasis + SWGBasis/PWCBasis/PWCHexBasis
# 全部在 P=1（单进程 MPI）下与串行装配逐项对照。
using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.Parallel
using LinearAlgebra
using MPI

if !MPI.Initialized()
    MPI.Init()
end
comm = MPI.COMM_WORLD

# ─────────────────────────────────────────────────────────────────────────────
# 小网格：单个立方体的四面体剖分 + 表面 + 单个 hex
# ─────────────────────────────────────────────────────────────────────────────
vol = generate_box_tet_mesh(0.4, 0.4, 0.4, 1, 1, 1)
surf = extract_surface(vol)
basis_surf = RWGBasis(surf)
basis_vol = SWGBasis(vol)
pwc = PWCBasis(vol)

hexnodes = Float64[
    0.0  1.0  1.0  0.0  0.0  1.0  1.0  0.0;
    0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0;
    0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
]
hexmesh = HexahedraMesh(1, hexnodes, reshape(collect(1:8), 8, 1), [1])
hexbasis = PWCHexBasis(hexmesh)

perms = fill(ComplexF64(2.0), num_elements(vol))
freq = 2e9
vefie = VEFIE(freq, perms)
scfie = SCFIE(freq, perms; alpha = 0.5)

function check_parallel_vs_serial(Z_mpi, Z_ser; rtol = 1e-6, label = "")
    Zf = gather(Z_mpi; root = 0)
    @test size(Zf) == size(Z_ser)
    @test norm(Zf - Z_ser) / norm(Z_ser) < rtol
    @test all(isfinite, Zf)
    if !isempty(label)
        @test norm(Zf) > 0
    end
end

@testset "VolumeAssembly: VEFIE + SWGBasis (P=1)" begin
    Z_mpi = assemble_impedance_matrix_parallel(vefie, basis_vol, perms)
    Z_ser = assemble_impedance_matrix(vefie, basis_vol)
    check_parallel_vs_serial(Z_mpi, Z_ser; label = "vefie-swg")
end

@testset "VolumeAssembly: VEFIE + PWCBasis (P=1)" begin
    # 显式 permittivities 路径
    Z_mpi = assemble_impedance_matrix_parallel(vefie, pwc, perms)
    Z_ser = assemble_impedance_matrix(vefie, pwc)
    check_parallel_vs_serial(Z_mpi, Z_ser; label = "vefie-pwc")
    # 默认 perms 包装路径
    Z_mpi2 = assemble_impedance_matrix_parallel(vefie, pwc)
    Zf2 = gather(Z_mpi2; root = 0)
    @test norm(Zf2 - Z_ser) / norm(Z_ser) < 1e-6
end

@testset "VolumeAssembly: VEFIE + PWCHexBasis (P=1)" begin
    Z_mpi = assemble_impedance_matrix_parallel(vefie, hexbasis)
    Z_ser = assemble_impedance_matrix(vefie, hexbasis)
    check_parallel_vs_serial(Z_mpi, Z_ser; label = "vefie-hex")
end

@testset "VolumeAssembly: SCFIE + RWGBasis + SWGBasis (P=1)" begin
    Z_mpi = assemble_impedance_matrix_parallel(scfie, basis_surf, basis_vol)
    Z_ser = assemble_impedance_matrix(scfie, basis_surf, basis_vol)
    check_parallel_vs_serial(Z_mpi, Z_ser; label = "scfie-swg")
end

@testset "VolumeAssembly: SCFIE + RWGBasis + PWCBasis (P=1)" begin
    Z_mpi = assemble_impedance_matrix_parallel(scfie, basis_surf, pwc)
    Z_ser = assemble_impedance_matrix(scfie, basis_surf, pwc)
    check_parallel_vs_serial(Z_mpi, Z_ser; label = "scfie-pwc")
end

@testset "VolumeAssembly: SCFIE + RWGBasis + PWCHexBasis (P=1)" begin
    Z_mpi = assemble_impedance_matrix_parallel(scfie, basis_surf, hexbasis)
    Z_ser = assemble_impedance_matrix(scfie, basis_surf, hexbasis)
    check_parallel_vs_serial(Z_mpi, Z_ser; label = "scfie-hex")
end
