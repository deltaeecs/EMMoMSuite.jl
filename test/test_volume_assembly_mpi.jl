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
        # VEFIE 对 lossless medium：实部 > 0（辐射阻抗）
        @test real(Z_full[1, 1]) > 0
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

end  # if isempty(mesh_file)
