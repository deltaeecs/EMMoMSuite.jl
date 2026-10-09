# test_sai_block.jl
# 逐八叉树块 SAI 预条件（SPAIPreconditioner(Z_near, cubes, row_map)）：
# 正确性 + 相比逐列 SAI 的迭代次数改善。合成块稀疏矩阵，fast 可跑。

using Test
using EMMoMSuite
using LinearAlgebra
using SparseArrays
using Random
using IterativeSolvers

using EMMoMSuite.Solvers: SPAIPreconditioner

# 仅需 bfInterval / neighbors 两个字段（Level.CubeInfo 的鸭子类型子集）
mutable struct TestCube
    bfInterval::UnitRange{Int}
    neighbors::Vector{Int}
end

@testset "SAI block preconditioner" begin
    Random.seed!(2026)

    # ── 1. 合成 4×4 八叉树叶层：16 cube，每 cube 4 个基函数 ──
    nSide = 4
    bfsPerCube = 4
    nCubes = nSide^2
    n = nCubes * bfsPerCube

    cubes = [TestCube(UnitRange((c-1)*bfsPerCube+1, c*bfsPerCube), Int[])
             for c = 1:nCubes]
    id3d = [(i, j) for i = 1:nSide for j = 1:nSide]
    for c = 1:nCubes
        i, j = id3d[c]
        for di = -1:1, dj = -1:1
            (di == 0 && dj == 0) && continue
            ii, jj = i+di, j+dj
            (1 <= ii <= nSide && 1 <= jj <= nSide) || continue
            push!(cubes[c].neighbors, (ii-1)*nSide + jj)
        end
        sort!(cubes[c].neighbors)
    end

    # ── 2. 按块近邻模式生成 Z_near（与 MLFMA 装配的 cube×neighbor 模式一致），
    #       对角占主但保留强近邻耦合 → 无预条件 GMRES 收敛慢 ──
    rows = Int[]; cols = Int[]; vals = ComplexF64[]
    for c = 1:nCubes
        neibfs = sort!(unique!(vcat([collect(cubes[i].bfInterval) for i in cubes[c].neighbors]...,
                                    collect(cubes[c].bfInterval))))
        for r in cubes[c].bfInterval, s in neibfs
            v = (randn() + 1im*randn()) * 0.5
            if r == s
                v = 10.0 + 2im   # 对角占主（EFIE 类对角强、近邻中强）
            end
            push!(rows, r); push!(cols, s); push!(vals, v)
        end
    end
    Z_near = sparse(rows, cols, vals, n, n)

    b = randn(ComplexF64, n)

    # ── 3. 构造：逐块 SAI（row_map = identity）与旧逐列 SAI ──
    P_block = SPAIPreconditioner(Z_near, cubes)               # 默认 row_map = identity
    P_col   = SPAIPreconditioner(Z_near)                      # 旧逐列路线（保持兼容）
    for (name, P) in [("block", P_block), ("column", P_col)]
        y = similar(b)
        ldiv!(y, P, b)
        @test all(isfinite, y)
        y2 = P \ b
        @test y2 ≈ y
    end

    # ── 4. 正确性：块 SAI 是比逐列更好的 A^{-1} 近似（M·Z ≈ I，左预条件方向）──
    Zd = Matrix(Z_near)
    I_eye = Matrix{ComplexF64}(I, n, n)
    res_block = opnorm(P_block.M * Zd - I_eye, 2)
    res_col   = opnorm(P_col.M * Zd   - I_eye, 2)
    @test res_block ≤ res_col * (1 + 1e-12)

    # ── 5. 迭代次数改善：GMRES(30) + 块 SAI 应少于无预条件，且不多于逐列 ──
    x0, ch0 = IterativeSolvers.gmres(Z_near, b; restart = 30, maxiter = 300,
                                     abstol = 1e-10, log = true)
    x1, ch1 = IterativeSolvers.gmres(Z_near, b; Pl = P_block, restart = 30,
                                     maxiter = 300, abstol = 1e-10, log = true)
    x2, ch2 = IterativeSolvers.gmres(Z_near, b; Pl = P_col, restart = 30,
                                     maxiter = 300, abstol = 1e-10, log = true)
    it0, it1, it2 = ch0.iters, ch1.iters, ch2.iters
    println("SAI block test — GMRES iters: none=$it0, column=$it2, block=$it1")
    @test it1 ≤ it0
    @test it1 ≤ it2
    @test norm(Z_near * x1 - b) / norm(b) ≤ 1e-8

    # ── 6. 非平凡 row_map（模拟 MLFMA sorted_ids 全局号索引）──
    perm = randperm(n)                       # sorted position -> global id
    # Z_perm[perm[i], perm[j]] == Z_near[i, j]，与 row_map = perm 一致
    ip = invperm(perm)
    Z_perm = Z_near[ip, ip]
    P_map = SPAIPreconditioner(Z_perm, cubes, i -> perm[i])
    y = similar(b)
    ldiv!(y, P_map, b)
    @test all(isfinite, y)
    # 置换一致：P_map.M 在全局号下应等于 P_block.M 在排序号下的结果
    @test P_map.M[perm, perm] ≈ P_block.M atol = 1e-8
end
