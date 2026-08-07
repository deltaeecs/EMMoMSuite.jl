"""
test_feko_reader.jl — Phase 14 TDD: Feko CSV 解析器单元测试

测试覆盖：
  14.0a  文件行结构 (1 header + 1442 data lines × 4 files)
  14.0b  θ/φ 范围与步长
  14.0c  RCS(m²) > 0 物理约束
  14.0d  RCS dBsm = 10*log10(sqm) 转换
  14.0e  split_phi_cuts 分组正确
  14.0f  文件不存在抛出合理错误
"""

using Test
using EMMoMSuite.Accuracy: read_feko_rcs, split_phi_cuts

const FEKO_DIR = joinpath(
    @__DIR__, "..", "deps", "fixtures", "AllinOne", "deps", "compare_feko"
)

const FEKO_FILES = [
    "jet_100MHzRCS.csv",
    "sphere_600MHzRCS.csv",
    "plate_1dot2GHzRCS.csv",
    "plate_metal_1dot2GHzRCS.csv",
]

@testset "14.0 Feko Reader" begin
    # ──────────────────────────────────────────────────────────────────────────
    # 14.0a  結構驗證 ─ 行數 & 各文件一致
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.0a 文件结构 (4 files, 1442 data lines each)" begin
        for fname in FEKO_FILES
            fpath = joinpath(FEKO_DIR, fname)
            if !isfile(fpath)
                @warn "Feko baseline not found, skipping: $fname"
                continue
            end
            theta, phi, rcs_sqm, rcs_dBsm = read_feko_rcs(fpath)
            @test length(theta) == 1442
            @test length(phi)   == 1442
            @test length(rcs_sqm) == 1442
            @test length(rcs_dBsm) == 1442
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.0b  θ/φ 范围与采样
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.0b θ/φ 范围" begin
        fpath = joinpath(FEKO_DIR, "jet_100MHzRCS.csv")
        isfile(fpath) || return
        theta, phi, _, _ = read_feko_rcs(fpath)

        # φ 只有两个切面
        phi_vals = sort(unique(round.(phi, digits = 1)))
        @test length(phi_vals) == 2
        @test phi_vals[1] ≈ 0.0   atol = 0.1
        @test phi_vals[2] ≈ 90.0  atol = 0.1

        # 每个切面 θ ∈ [-180, 180]
        for phi_cut in phi_vals
            mask = abs.(phi .- phi_cut) .< 0.1
            th = theta[mask]
            @test minimum(th) ≈ -180.0  atol = 0.1
            @test maximum(th) ≈  180.0  atol = 0.1
            @test length(th) == 721
        end

        # 步长 0.5°
        th0 = sort(theta[abs.(phi) .< 0.1])
        diffs = diff(th0)
        @test all(abs.(diffs .- 0.5) .< 1e-6)
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.0c  RCS 物理约束
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.0c RCS > 0 (物理约束)" begin
        for fname in FEKO_FILES
            fpath = joinpath(FEKO_DIR, fname)
            isfile(fpath) || continue
            _, _, rcs_sqm, _ = read_feko_rcs(fpath)
            @test all(rcs_sqm .> 0.0)
        end
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.0d  dBsm 转换正确
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.0d dBsm = 10*log10(sqm)" begin
        fpath = joinpath(FEKO_DIR, "jet_100MHzRCS.csv")
        isfile(fpath) || return
        _, _, rcs_sqm, rcs_dBsm = read_feko_rcs(fpath)
        expected = 10.0 .* log10.(rcs_sqm)
        @test maximum(abs.(rcs_dBsm .- expected)) < 1e-10
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.0e  split_phi_cuts 分组
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.0e split_phi_cuts" begin
        fpath = joinpath(FEKO_DIR, "sphere_600MHzRCS.csv")
        isfile(fpath) || return
        theta, phi, _, rcs_dBsm = read_feko_rcs(fpath)
        cuts = split_phi_cuts(theta, phi, rcs_dBsm)

        @test haskey(cuts, 0.0)
        @test haskey(cuts, 90.0)
        @test length(cuts[0.0].theta)    == 721
        @test length(cuts[90.0].theta)   == 721
        @test length(cuts[0.0].rcs_dBsm) == 721

        # 两个切面 θ 有序
        @test issorted(cuts[0.0].theta)
        @test issorted(cuts[90.0].theta)
    end

    # ──────────────────────────────────────────────────────────────────────────
    # 14.0f  文件不存在时抛出 ArgumentError
    # ──────────────────────────────────────────────────────────────────────────
    @testset "14.0f 错误处理" begin
        @test_throws ArgumentError read_feko_rcs("/nonexistent/path.csv")
    end
end
