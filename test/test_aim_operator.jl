# AIM/IE-FFT（P3）测试：FFT 卷积正确性、投影一致性、matvec vs 稠密参照
using Test
using EMMoMSuite
using LinearAlgebra, Random, SparseArrays
using StaticArrays
import EMMoMSuite.FastAlgorithms.AIM as AIM
using EMMoMSuite.Geometry, EMMoMSuite.BasisFunctions, EMMoMSuite.IntegralEquations
using EMMoMSuite.IntegralEquations.EFIEModule: efie_interaction!

@testset "AIM/IE-FFT (P3)" begin

    @testset "FFT 卷积 = 直接求和（机器精度）" begin
        h = 0.1
        g = AIM.AIMGrid{Float64}(SVector(0.0, 0.0, 0.0), (6, 6, 6), h)
        k = 2π
        Ghat = AIM.build_kernel(g, k)
        P = size(Ghat)
        n = g.n
        NP = prod(n)
        Random.seed!(5)
        U = randn(ComplexF64, NP)
        V = zeros(ComplexF64, NP)
        AIM.conv3!(V, U, Ghat, P, n)
        # 直接求和
        Vd = zeros(ComplexF64, NP)
        for i in 1:NP, j in 1:NP
            Vd[i] += AIM.greens(k, norm(AIM.node_coord(g, i) - AIM.node_coord(g, j))) * U[j]
        end
        @test norm(V - Vd) / norm(Vd) < 1e-10
    end

    @testset "投影一致性（零阶矩与散度）" begin
        freq = 300e6
        λ = 299792458.0 / freq
        mesh = generate_sphere_mesh(0.2, 6, 12)
        basis = RWGBasis(mesh)
        g = AIM.grid_from_mesh(mesh, 0.1 * λ)
        tris = EMMoMSuite.IntegralEquations.Impedance.get_triangles_info(mesh, basis)
        nodes, Pv, Pd = AIM.build_projection(basis, g, tris)
        for m in eachindex(basis.functions)
            # Σ Πv = ∫ f dS（三线性 Λ 完备性 ΣΛ=1）
            intf = sum(Pv[m])
            bf = basis.functions[m]
            ref = zeros(3)
            for i_supp in 1:2
                tri_idx = bf.support[i_supp]
                tri_idx == 0 && continue
                verts = tris[tri_idx].vertices
                v_opp = verts[:, bf.local_edge_idx[i_supp]]
                for q in 1:7
                    r = verts * AIM._TRI7_L[q]
                    ref .+= AIM._TRI7_W[q] * (bf.signs[i_supp] * bf.edge_length / 2) .* (r - v_opp)
                end
            end
            @test norm(intf - ref) / max(norm(ref), 1e-30) < 1e-10
            # Σ Πd = ∫ ∇·f dS = 0（RWG 无散度净通量）
            @test abs(sum(Pd[m])) < 1e-10
        end
    end

    @testset "matvec vs 稠密 Z（EFIE，0.8λ 球，走远场 FFT 路径）" begin
        freq = 300e6
        λ = 299792458.0 / freq
        mesh = generate_sphere_mesh(0.8, 14, 28)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        N = num_basis(basis)
        Z_dense = assemble_impedance_matrix(efie, basis)
        op = AIMOperator(efie, basis; h_ratio = 0.1)
        @test op isa AIMOperator
        @test size(op) == (N, N)
        Random.seed!(8)
        x = randn(ComplexF64, N)
        y_aim = op * x
        y_dense = Z_dense * x
        rel = norm(y_aim - y_dense) / norm(y_dense)
        @info "AIM vs dense matvec rel err" rel N grid = op.grid.n
        @test rel < 1e-2
        # 近场修正必须真实存在（非全部精确对）
        @test nnz(op.Z_near) < N * N
    end

    @testset "近场装配空间哈希：候选覆盖 = O(N²) 重叠模式 + 值抽查" begin
        freq = 300e6
        λ = 299792458.0 / freq
        mesh = generate_sphere_mesh(0.3, 8, 16)
        basis = RWGBasis(mesh)
        efie = EFIE(freq)
        op = AIMOperator(efie, basis; h_ratio = 0.1, near_radius = 0.35)
        N = num_basis(basis)
        tris = EMMoMSuite.IntegralEquations.Impedance.get_triangles_info(mesh, basis)
        k = op.k
        nr = 0.35 * λ
        boxes = Vector{NTuple{6,Float64}}(undef, N)
        for m in 1:N
            bf = basis.functions[m]
            lo = fill(Inf, 3); hi = fill(-Inf, 3)
            for s in 1:2
                t = bf.support[s]
                t == 0 && continue
                for v in eachcol(tris[t].vertices), d in 1:3
                    lo[d] = min(lo[d], v[d]); hi[d] = max(hi[d], v[d])
                end
            end
            boxes[m] = (lo[1]-nr, hi[1]+nr, lo[2]-nr, hi[2]+nr, lo[3]-nr, hi[3]+nr)
        end
        overlap(b1, b2) =
            b1[1] <= b2[2] && b2[1] <= b1[2] &&
            b1[3] <= b2[4] && b2[3] <= b1[4] &&
            b1[5] <= b2[6] && b2[5] <= b1[6]
        # 1) 非零模式：哈希装配的 Z_near 非零模式必须与 O(N²) 盒重叠模式一致（含自对）
        Zsp = sparse(op.Z_near)
        for m in 1:N, n in 1:N
            expect = m == n || overlap(boxes[m], boxes[n])
            @test (Zsp[m, n] != 0) == expect
        end
        # 2) 值抽查：对角线 + 一个近邻对，与 Z_direct − Z_grid 一致
        C = im * k * op.eta
        for (m, n) in ((1, 1), (2, 1), (1, 3))
            bf_m = basis.functions[m]; bf_n = basis.functions[n]
            Zdir = zero(ComplexF64)
            Z3 = zeros(ComplexF64, 3, 3)
            for a in 1:2, b in 1:2
                tm = bf_m.support[a]; sn = bf_n.support[b]
                (tm == 0 || sn == 0) && continue
                fill!(Z3, 0)
                efie_interaction!(Z3, efie, tris[tm], tris[sn])
                Zdir +=
                    Z3[bf_m.local_edge_idx[a], bf_n.local_edge_idx[b]] *
                    bf_m.signs[a] * bf_n.signs[b]
            end
            Zgrid = zero(ComplexF64)
            for (qi, ni) in enumerate(op.nodes[m]), (qj, nj) in enumerate(op.nodes[n])
                G = AIM.greens(k, norm(AIM.node_coord(op.grid, ni) - AIM.node_coord(op.grid, nj)))
                Zgrid += dot(op.Pv[m][qi], op.Pv[n][qj]) * G
                Zgrid -= op.Pd[m][qi] * op.Pd[n][qj] * G / k^2
            end
            @test abs(op.Z_near[m, n] - (Zdir - C * Zgrid)) < 1e-10
        end
    end

end
