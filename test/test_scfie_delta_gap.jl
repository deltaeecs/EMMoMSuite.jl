"""
test_scfie_delta_gap.jl — Phase 15 步骤 15.4 / 15.5

测试覆盖：
  15.4/15.5  excitation_vector(op::SCFIE, source::DeltaGapSource,
                               surf_basis::RWGBasis, vol_basis::SWGBasis)
             → 长度 (N_surf + N_vol) 向量；
               表面部分 (1:N_surf) 有 delta-gap；体积部分全零
"""

using Test
using EMMoMSuite
using LinearAlgebra

@testset "15.4/15.5 SCFIE DeltaGap 激励" begin

    # ── 查找 TriTetra.nas（混合表面+体积网格） ─────────────────────────
    mesh_file = ""
    for candidate in [
        joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne", "meshfiles", "TriTetra.nas"),
        joinpath(@__DIR__, "..", "deps", "fixtures", "Basics",   "meshfiles", "TriTetra.nas"),
        joinpath(@__DIR__, "..", "deps", "fixtures", "Kernels",  "meshfiles", "TriTetra.nas"),
    ]
        if isfile(candidate)
            mesh_file = candidate
            break
        end
    end

    if isempty(mesh_file)
        @warn "TriTetra.nas 未找到 — 跳过 SCFIE DeltaGap 测试（步骤 15.4/15.5）"
        @test_skip true
    else
        surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale = 0.001)
        surf_basis = RWGBasis(surf_mesh)
        vol_basis  = SWGBasis(vol_mesh)
        N_surf = num_basis(surf_basis)
        N_vol  = num_basis(vol_basis)

        @test N_surf > 0
        @test N_vol  > 0

        freq  = 1e9
        eps_r = 2.0
        perms = fill(ComplexF64(eps_r), num_elements(vol_mesh))
        scfie = SCFIE(freq, perms; alpha = 0.5)

        V_feed = 1.0 + 0im
        feed   = DeltaGapSource(freq, [1], V_feed)   # 第 1 条表面 RWG 边
        ℓ₁     = surf_basis.functions[1].edge_length

        # ─────────────────────────────────────────────────────────────────
        # 15.5: excitation_vector(SCFIE, DeltaGapSource, RWGBasis, SWGBasis)
        # ─────────────────────────────────────────────────────────────────
        @testset "vector length is N_surf + N_vol" begin
            V = excitation_vector(scfie, feed, surf_basis, vol_basis)
            @test length(V) == N_surf + N_vol
        end

        @testset "surface feed edge has correct value" begin
            V = excitation_vector(scfie, feed, surf_basis, vol_basis)
            @test V[1] ≈ V_feed * ℓ₁
        end

        @testset "surface non-feed entries are zero" begin
            V = excitation_vector(scfie, feed, surf_basis, vol_basis)
            if N_surf > 1
                @test all(iszero, V[2:N_surf])
            end
        end

        @testset "volume part is all zeros" begin
            V = excitation_vector(scfie, feed, surf_basis, vol_basis)
            V_vol = V[N_surf+1:end]
            @test all(iszero, V_vol)
        end

        @testset "multiple surface feed edges" begin
            if N_surf >= 2
                feed2 = DeltaGapSource(freq, [1, 2], V_feed)
                V     = excitation_vector(scfie, feed2, surf_basis, vol_basis)
                ℓ₂    = surf_basis.functions[2].edge_length
                @test V[1] ≈ V_feed * ℓ₁
                @test V[2] ≈ V_feed * ℓ₂
                @test all(iszero, V[N_surf+1:end])
            end
        end

        @testset "out-of-range surface edge is skipped" begin
            bad  = DeltaGapSource(freq, [N_surf + 1], V_feed)   # 超出表面基函数范围
            V    = excitation_vector(scfie, bad, surf_basis, vol_basis)
            @test all(iszero, V)
        end
    end
end
