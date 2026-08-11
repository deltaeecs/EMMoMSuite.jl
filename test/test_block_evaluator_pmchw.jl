# test/test_block_evaluator_pmchw.jl - F2 PMCHW 2N block evaluator (TDD)
#
# PMCHWBlockEvaluator 按 2N 系统全局索引求值远场块 Z[rows, cols]（行/列可混合
# J/M 通道）。生产用法：ACA/MLACA 只对八叉树非邻盒子对（远场）调用本求值器，
# 近场（自对/邻对，含半解析近场项的顺序不对称性）由 Z_near 精确装配覆盖。
# 测试取 octree 非邻叶子对的 J/M 索引与稠密 PMCHW 子矩阵逐元素对照。

using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.FastAlgorithms.ACA
using EMMoMSuite.FastAlgorithms.MLFMA
using LinearAlgebra
using Random

@testset "PMCHW block evaluator" begin
    mesh = generate_sphere_mesh(0.5, 6, 10)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    pmchw = PMCHW(300e6, 4.0)
    Z = assemble_impedance_matrix(pmchw, basis)

    ev = PMCHWBlockEvaluator(pmchw, basis)

    # 取八叉树两个非邻叶子盒子的 J/M 全局索引作为远场块
    bf_centers = reduce(hcat, [bf.center for bf in basis.functions])
    λ = 299792458.0 / 300e6
    octree, sorted_ids = build_octree(bf_centers, 0.25 * λ; λ = λ, near_range = 1)
    leaf = octree.levels[octree.nLevels]
    cubes = leaf.cubes
    pair = nothing
    for i in eachindex(cubes)
        isempty(cubes[i].bfInterval) && continue
        for j in (i + 1):length(cubes)
            isempty(cubes[j].bfInterval) && continue
            j in cubes[i].neighbors && continue
            pair = (i, j)
            break
        end
        pair !== nothing && break
    end
    @test pair !== nothing
    i, j = pair
    rowsJ = sorted_ids[cubes[i].bfInterval]
    colsJ = sorted_ids[cubes[j].bfInterval]
    rows = vcat(rowsJ, N .+ rowsJ)
    cols = vcat(colsJ, N .+ colsJ)

    B = eval_block(ev, rows, cols)
    @test size(B) == (length(rows), length(cols))
    for (ii, i) in enumerate(rows), (jj, j) in enumerate(cols)
        @test isapprox(B[ii, jj], Z[i, j]; atol = 1e-10)
    end

    # 更大的远场块：从两个非邻盒子的并集中随机取混合 J/M 行/列
    rng = Random.MersenneTwister(3)
    pool_r = vcat(rows, rows .+ 0)  # 保持 J/M 混合
    rows2 = sort(unique(rand(rng, pool_r, 24)))
    cols2 = sort(unique(rand(rng, vcat(cols, cols .+ 0), 24)))
    B2 = eval_block(ev, rows2, cols2)
    @test size(B2) == (length(rows2), length(cols2))
    for (ii, i) in enumerate(rows2), (jj, j) in enumerate(cols2)
        @test isapprox(B2[ii, jj], Z[i, j]; atol = 1e-10)
    end
end
