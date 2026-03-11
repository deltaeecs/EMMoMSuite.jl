# test_volume_assembly_mpi.jl
# 覆盖率目标: VolumeAssembly.jl — VEFIE+SWG 和 SCFIE+RWG+SWG 并行装配
# 使用 P=1 单进程直接调用，无需 mpiexec，覆盖率工具可正常捕获

using Test
using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.Parallel
using LinearAlgebra
using MPI

# TriTetra.nas 路径解析（与 test_scfie.jl 相同的查找逻辑）
function find_tritetra_mesh()
    candidates = [
        joinpath(@__DIR__, "..", "..", "MoM_Basics",  "meshfiles", "TriTetra.nas"),
        joinpath(@__DIR__, "..", "..", "MoM_Kernels", "meshfiles", "TriTetra.nas"),
    ]
    for c in candidates
        isfile(c) && return c
    end
    return ""
end

mesh_file = find_tritetra_mesh()

if isempty(mesh_file)
    @warn "TriTetra.nas not found — skipping VolumeAssembly MPI tests"
    @testset "VolumeAssembly MPI (skipped)" begin @test_skip true end
else

# ─────────────────────────────────────────────────────────────────────────────
# 共用设置
# ─────────────────────────────────────────────────────────────────────────────
if !MPI.Initialized()
    MPI.Init()
end
init_parallel!()

surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001)
basis_surf = RWGBasis(surf_mesh)
basis_vol  = SWGBasis(vol_mesh)

freq   = 2e9
eps_r  = complex(2.0, -0.001)
perms  = fill(eps_r, num_elements(vol_mesh))
scfie  = SCFIE(freq, perms; alpha=0.5)
vefie  = VEFIE(freq, perms)

n_surf  = num_basis(basis_surf)
n_vol   = num_basis(basis_vol)
n_total = n_surf + n_vol

# ─────────────────────────────────────────────────────────────────────────────
# VEFIE + SWGBasis 并行装配
# ─────────────────────────────────────────────────────────────────────────────
@testset "VolumeAssembly: VEFIE + SWGBasis (P=1)" begin
    Z_mpi = assemble_impedance_matrix_parallel(vefie, basis_vol, perms)

    @test size(Z_mpi) == (n_vol, n_vol)

    Z_full = gather(Z_mpi; root=0)
    if Z_full !== nothing
        # Z[1,1] 应为有限复数
        @test isfinite(real(Z_full[1, 1]))
        @test isfinite(imag(Z_full[1, 1]))
        # VEFIE 对 lossless medium：实部应接近非负，允许零附近的数值摆动
        @test real(Z_full[1, 1]) > -1e-6
        # 非零矩阵
        @test norm(Z_full) > 0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# SCFIE + RWGBasis + SWGBasis 并行装配
# ─────────────────────────────────────────────────────────────────────────────
@testset "VolumeAssembly: SCFIE + RWGBasis + SWGBasis (P=1)" begin
    Z_mpi = assemble_impedance_matrix_parallel(scfie, basis_surf, basis_vol)

    @test size(Z_mpi) == (n_total, n_total)

    Z_full = gather(Z_mpi; root=0)
    if Z_full !== nothing
        # 对比已知基准值 (test_mpi_scfie_small.jl P=1 输出)
        @test isapprox(Z_full[1, 1], 0.010743309358896107 - 0.03340102794592738im; rtol=1e-5)
        @test norm(Z_full) > 0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# VEFIE + PWCBasis 并行装配
# ─────────────────────────────────────────────────────────────────────────────
@testset "VolumeAssembly: VEFIE + PWCBasis (P=1)" begin
    basis_pwc = PWCBasis(vol_mesh)
    n_pwc = num_basis(basis_pwc)   # 3 × n_tets

    Z_mpi = assemble_impedance_matrix_parallel(vefie, basis_pwc)

    @test size(Z_mpi) == (n_pwc, n_pwc)

    Z_full = gather(Z_mpi; root=0)
    if Z_full !== nothing
        @test isfinite(real(Z_full[1, 1]))
        @test real(Z_full[1, 1]) > 0
        @test norm(Z_full) > 0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# SCFIE + RWGBasis + PWCBasis 并行装配
# ─────────────────────────────────────────────────────────────────────────────
@testset "VolumeAssembly: SCFIE + RWGBasis + PWCBasis (P=1)" begin
    basis_pwc  = PWCBasis(vol_mesh)
    n_pwc      = num_basis(basis_pwc)
    n_tot_pwc  = n_surf + n_pwc

    Z_mpi = assemble_impedance_matrix_parallel(scfie, basis_surf, basis_pwc)

    @test size(Z_mpi) == (n_tot_pwc, n_tot_pwc)

    Z_full = gather(Z_mpi; root=0)
    if Z_full !== nothing
        @test isfinite(real(Z_full[1, 1]))
        @test norm(Z_full) > 0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# VEFIE + PWCHexBasis 并行装配（内联构造六面体网格，无需 .nas 文件）
# ─────────────────────────────────────────────────────────────────────────────
@testset "VolumeAssembly: VEFIE + PWCHexBasis (P=1)" begin
    # 构造 2 个相邻六面体 cube: x∈[0,1] 和 x∈[1,2]，y/z∈[0,1]
    _nodes_hex = Float64[
        0.0  1.0  1.0  0.0  0.0  1.0  1.0  0.0   2.0  2.0  2.0  2.0;
        0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0   0.0  1.0  1.0  0.0;
        0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0   0.0  0.0  1.0  1.0
    ]
    _els_hex = Int[1 2; 2 9; 3 10; 4 3; 5 6; 6 11; 7 12; 8 7]
    _hex_mesh = HexahedraMesh(2, _nodes_hex, _els_hex, [1, 1])

    _perms_hex = fill(complex(2.0, -0.001), 2)   # 2 hexahedra
    _vefie_hex = VEFIE(2e9, _perms_hex)
    _basis_pwchex = PWCHexBasis(_hex_mesh)
    _n_pwchex = num_basis(_basis_pwchex)

    Z_mpi_hex = assemble_impedance_matrix_parallel(_vefie_hex, _basis_pwchex)

    @test size(Z_mpi_hex) == (_n_pwchex, _n_pwchex)

    Z_full_hex = gather(Z_mpi_hex; root=0)
    if Z_full_hex !== nothing
        @test isfinite(real(Z_full_hex[1, 1]))
        @test isfinite(imag(Z_full_hex[1, 1]))
        @test norm(Z_full_hex) > 0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# SCFIE + RWGBasis + PWCHexBasis 并行装配（内联最小网格）
# ─────────────────────────────────────────────────────────────────────────────
@testset "VolumeAssembly: SCFIE + RWGBasis + PWCHexBasis (P=1)" begin
    # 表面网格：2个三角形，1条公共内边 → 1个 RWG DOF
    #   节点（小尺寸 mm 量级，避免 numerical overflow）
    _nodes_tri = Float64[
        0.0   1.0   0.0   1.0;
        0.0   0.0   1.0   1.0;
        0.0   0.0   0.0   0.0
    ] * 1e-3   # 1mm
    _els_tri = [1 2; 2 4; 3 3]   # T1=1-2-3, T2=2-4-3; 共享边 2-3
    _tri_mesh = TriangleMesh(2, _nodes_tri, _els_tri, [1, 1])
    _basis_rwg = RWGBasis(_tri_mesh)

    # 六面体网格：单个 1mm³ 六面体
    _scale = 1e-3
    _nodes_hex2 = Float64[
        0.0  1.0  1.0  0.0  0.0  1.0  1.0  0.0;
        0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0;
        0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
    ] * _scale
    _els_hex2 = reshape(collect(1:8), 8, 1)
    _hex_mesh2 = HexahedraMesh(1, _nodes_hex2, _els_hex2, [1])

    _perms_hex2   = fill(complex(2.0, -0.001), 1)   # 1 hexahedra
    _scfie_hex    = SCFIE(2e9, _perms_hex2; alpha=0.5)
    _basis_pwchex2 = PWCHexBasis(_hex_mesh2)
    _n_pwchex2    = num_basis(_basis_pwchex2)   # 3
    _n_rwg        = num_basis(_basis_rwg)         # 1
    _n_tot_h      = _n_rwg + _n_pwchex2

    Z_mpi_sh = assemble_impedance_matrix_parallel(_scfie_hex, _basis_rwg, _basis_pwchex2)

    @test size(Z_mpi_sh) == (_n_tot_h, _n_tot_h)

    Z_full_sh = gather(Z_mpi_sh; root=0)
    if Z_full_sh !== nothing
        @test isfinite(real(Z_full_sh[1, 1]))
        @test norm(Z_full_sh) > 0
    end
end

end  # if isempty(mesh_file)
