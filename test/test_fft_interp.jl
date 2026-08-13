# FFT φ 方向谱插值（方案 P1）测试
# 门：带限函数 FFT 插值达机器精度；与 Lagrange ϕCSC 结果一致；反插值往返成立；MLFMA 集成不劣化。

using Test
using EMMoMSuite
using LinearAlgebra, Random, SparseArrays
import EMMoMSuite.FastAlgorithms.MLFMA.Interpolation as Interp
using EMMoMSuite.FastAlgorithms.MLFMA: MLFMAOperator
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations

@testset "FFT φ 谱插值 (P1)" begin

    @testset "带限函数 FFT 插值 = 精确采样（机器精度）" begin
        M1, M2, nθ = 32, 64, 5
        K = 10   # 带宽 K < M1/2，无混叠
        Random.seed!(7)
        coeffs = randn(ComplexF64, nθ, 2K + 1)
        x = zeros(ComplexF64, nθ * M1)
        for j in 1:M1, i in 1:nθ
            φ = 2π * (j - 0.5) / M1   # MLFMA 半格偏置网格
            x[(j-1)*nθ + i] = sum(coeffs[i, m+K+1] * exp(im * m * φ) for m in -K:K)
        end
        y = Interp.fft_interp_phi(x, nθ, M1, M2)
        err = 0.0
        for j in 1:M2, i in 1:nθ
            φ = 2π * (j - 0.5) / M2
            ex = sum(coeffs[i, m+K+1] * exp(im * m * φ) for m in -K:K)
            err = max(err, abs(y[(j-1)*nθ + i] - ex))
        end
        @test err < 1e-9
    end

    @testset "与 Lagrange ϕCSC 结果一致（光滑函数）" begin
        Lc, poles_c = Interp.levelIntegralInfoCal(1.0; λ = 1.0, L_min = 0)
        Lp, poles_p = Interp.levelIntegralInfoCal(2.0; λ = 1.0, L_min = 0)
        nθ = length(poles_c.Xθs)
        M1 = length(poles_c.Xϕs)
        M2 = length(poles_p.Xϕs)
        @test M2 > M1

        f(θ, φ) = exp(0.3 * cos(φ)) * (2 + cos(θ))
        x = zeros(ComplexF64, nθ * M1)
        for j in 1:M1, i in 1:nθ
            x[(j-1)*nθ + i] = f(poles_c.Xθs[i], poles_c.Xϕs[j])   # (θ 内、φ 外)，φ 为半格偏置
        end

        info = Interp.interpolationCSCMatCal(poles_p, poles_c, 6)
        y_lag = info.ϕCSC * x
        y_fft = Interp.fft_interp_phi(x, nθ, M1, M2)

        @test length(y_lag) == nθ * M2
        @test length(y_fft) == nθ * M2
        rel = maximum(abs.(y_fft - y_lag)) / maximum(abs.(y_lag))
        @test rel < 1e-4
        # 与精确值对比（φ 方向 exp(0.3cosφ) 带宽远小于 Nyquist）
        exact = zeros(ComplexF64, nθ * M2)
        for j in 1:M2, i in 1:nθ
            exact[(j-1)*nθ + i] = f(poles_c.Xθs[i], poles_p.Xϕs[j])
        end
        @test maximum(abs.(y_fft - exact)) < 1e-8
        @test maximum(abs.(y_lag - exact)) < 1e-3
    end

    @testset "反插值 = 插值矩阵的转置（与 Lagrange ϕCSCT 语义一致）" begin
        M1, M2 = 8, 16
        Random.seed!(11)
        # 单 θ 行：显式构造插值矩阵 P 与反插值矩阵 Q，验证 Q == P'
        P = zeros(ComplexF64, M2, M1)
        for j in 1:M1
            e = zeros(ComplexF64, M1); e[j] = 1.0
            P[:, j] = Interp.fft_interp_phi(e, 1, M1, M2)
        end
        Q = zeros(ComplexF64, M1, M2)
        for j in 1:M2
            e = zeros(ComplexF64, M2); e[j] = 1.0
            Q[:, j] = Interp.fft_anterp_phi(e, 1, M2, M1)
        end
        @test norm(Q - transpose(P)) / norm(P) < 1e-12
        # 多 θ 行同样成立（分块对角结构）
        nθ = 3
        P3 = zeros(ComplexF64, nθ * M2, nθ * M1)
        for j in 1:(nθ*M1)
            e = zeros(ComplexF64, nθ * M1); e[j] = 1.0
            P3[:, j] = Interp.fft_interp_phi(e, nθ, M1, M2)
        end
        Q3 = zeros(ComplexF64, nθ * M1, nθ * M2)
        for j in 1:(nθ*M2)
            e = zeros(ComplexF64, nθ * M2); e[j] = 1.0
            Q3[:, j] = Interp.fft_anterp_phi(e, nθ, M2, M1)
        end
        @test norm(Q3 - transpose(P3)) / norm(P3) < 1e-12
    end

    @testset "构造 FFTSpectral 极点信息与插值类型" begin
        L, poles = Interp.levelIntegralInfoCal(1.0, Val(:FFTSpectral); λ = 1.0, L_min = 0)
        @test poles isa Interp.FFTGLPolesInfo
        @test Interp.interp_type(poles) === Interp.FFTInterpInfo{Int,Float64}
        @test length(poles.Xθs) > 0 && length(poles.Xϕs) > 0
    end

    @testset "批量 FFT 插值/反插值 = 逐调用结果" begin
        Lc, poles_c = Interp.levelIntegralInfoCal(0.5; λ = 1.0, L_min = 0)
        Lp, poles_p = Interp.levelIntegralInfoCal(1.0; λ = 1.0, L_min = 0)
        info_lag = Interp.interpolationCSCMatCal(poles_p, poles_c, 6)
        info = Interp.FFTInterpInfo(
            info_lag.θCSC,
            info_lag.θCSCT,
            length(poles_c.Xθs),
            length(poles_c.Xϕs),
            length(poles_p.Xϕs),
            7,
        )
        nCubes = 7
        Random.seed!(21)
        agg = randn(ComplexF64, info.nθ * info.M1, 2, nCubes)
        out = zeros(ComplexF64, info.nθ * info.M2, 2, nCubes)
        Interp.fft_interp_phi_batch!(out, agg, info)
        for c in 1:nCubes, p in 1:2
            ref = Interp.fft_interp_phi(vec(agg[:, p, c]), info.nθ, info.M1, info.M2)
            @test maximum(abs.(out[:, p, c] .- ref)) < 1e-10
        end
        temp = randn(ComplexF64, info.nθ * info.M2, 2, nCubes)
        out2 = zeros(ComplexF64, info.nθ * info.M1, 2, nCubes)
        Interp.fft_anterp_phi_batch!(out2, temp, info)
        for c in 1:nCubes, p in 1:2
            ref = Interp.fft_anterp_phi(vec(temp[:, p, c]), info.nθ, info.M2, info.M1)
            @test maximum(abs.(out2[:, p, c] .- ref)) < 1e-10
        end
    end

    @testset "MLFMA 集成：FFTSpectral 与 Lagrange2Step matvec 一致" begin
        freq = 300e6
        λ = 299792458.0 / freq
        # 半径 1.0λ 的大球 + 0.2λ 叶层 → 4 层，粗层存在远亲盒，远场被实际执行
        mesh = generate_sphere_mesh(1.0, 12, 24)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        leaf = 0.2 * λ

        op_lag = MLFMAOperator(efie, basis, leaf, Val(:Lagrange2Step))
        op_fft = MLFMAOperator(efie, basis, leaf, Val(:FFTSpectral))

        @test op_lag.octree.nLevels >= 4
        @test op_fft.octree.nLevels >= 4
        @test op_lag.octree.levels[op_lag.octree.nLevels-1].interpWθϕ isa EMMoMSuite.FastAlgorithms.MLFMA.Interpolation.LagrangeInterpInfo
        @test op_fft.octree.levels[op_fft.octree.nLevels-1].interpWθϕ isa Interp.FFTInterpInfo

        Random.seed!(3)
        x = randn(ComplexF64, num_basis(basis))
        Z_dense = assemble_impedance_matrix(efie, basis)
        y_lag = op_lag * x
        y_fft = op_fft * x
        y_dense = Z_dense * x

        # 远场必须实际有贡献，否则插值差异不会体现在 matvec 上
        @test norm(y_lag - op_lag.Z_near * x) / norm(y_lag) > 1e-6
        # FFT 谱插值相对稠密参照的误差不得劣于 Lagrange 2 倍以上
        err_lag = norm(y_lag - y_dense) / norm(y_dense)
        err_fft = norm(y_fft - y_dense) / norm(y_dense)
        @test err_fft < 2 * err_lag
        @info "FFTSpectral vs Lagrange2Step vs dense" err_fft err_lag
        rel = norm(y_fft - y_lag) / norm(y_lag)
        @test rel < 1e-2
        @info "FFTSpectral vs Lagrange2Step matvec rel diff" rel
    end

end
