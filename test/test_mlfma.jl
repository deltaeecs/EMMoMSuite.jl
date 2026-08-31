using Test
using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.CoreModule
using EMMoMSuite.Solvers: BlockJacobiPreconditioner
using StaticArrays
using LinearAlgebra
using SpecialFunctions
using MPI

# Mock Operator
struct MockOperator <: AbstractIntegralOperator
    k::Float64
end

struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

@testset "MLFMA System" begin
    # 1. Setup Geometry and Basis
    # Create a larger mesh (2x2 grid of squares, 8 triangles) to ensure enough basis functions
    # Vertices
    x = [0.0, 1.0, 2.0]
    y = [0.0, 1.0, 2.0]
    nodes = zeros(3, 9)
    k = 1
    for j in 1:3, i in 1:3
        nodes[1, k] = x[i]
        nodes[2, k] = y[j]
        k += 1
    end
    
    # Elements
    elements = Int[]
    # Node indices in grid:
    # 7 8 9
    # 4 5 6
    # 1 2 3
    # (Note: my loop above generates 1,2,3 for y=0, etc. which matches this visual if y increases upwards)
    
    for j in 1:2, i in 1:2
        n1 = (j-1)*3 + i
        n2 = n1 + 1
        n3 = n1 + 3
        n4 = n1 + 4
        # Quad n1-n2-n4-n3
        # Split into two triangles: n1-n2-n4 and n1-n4-n3
        push!(elements, n1, n2, n4)
        push!(elements, n1, n4, n3)
    end
    elements = reshape(elements, 3, :)
    
    tags = ones(Int, size(elements, 2))
    mesh = TriangleMesh(size(elements, 2), nodes, elements, tags)
    basis = RWGBasis(mesh)
    
    # Get basis centers for octree construction
    # For RWG, center is usually edge center.
    # We must use basis function centers so that the octree sorts the basis functions.
    bf_centers = reduce(hcat, [bf.center for bf in basis.functions])
    
    leafCubeEdgel = 0.1 # Small enough to create multiple levels
    
    # 2. Build Octree
    octree, sorted_ids = build_octree(bf_centers, leafCubeEdgel; λ=1.0)
    
    # Reorder basis functions to match the octree sorting
    permute!(basis.functions, sorted_ids)
    
    @test octree isa OctreeInfo
    @test octree.nLevels >= 2
    
    # 3. Precomputations
    # Check if precomputations ran during build_octree
    # Shift factors
    level2 = octree.levels[2]
    if octree.nLevels > 2
        @test !isempty(level2.phaseShift2Kids)
    end
    
    # Interpolation matrices
    # Should be present in child levels (e.g. leaf level if nLevels > 2)
    # Note: build_octree calls compute_interpolation_matrices!
    
    # Translation factors
    # build_octree calls compute_translation_factors!；αTrans 只保存实际远邻偏移列，
    # 无 farneighbors 的层（小八叉树层级 2）允许为空表（translate! 语义不变）。
    @test isdefined(level2, :αTrans)
    @test isdefined(level2, :αTransIndex)
    
    # 4. Aggregation
    operator = MockOperator(1.0) # k=1.0
    
    # Run aggregation
    x = ones(ComplexF64, length(basis.functions))
    aggregate!(octree, basis, operator, x, sorted_ids)
    
    # Check if aggregation arrays are populated
    leafLevel = octree.levels[octree.nLevels]
    @test !isempty(leafLevel.aggS)
    # Check if values are not all zero (assuming basis functions exist and k!=0)
    # With 1 basis function, it should have some value.
    # However, aggregate_leaf! logic depends on basis function support being in the cube.
    # We need to ensure the basis function is associated with a cube.
    # The octree was built on vertices. The basis function centers might not exactly match vertices, 
    # but setBFInterval! (called in build_octree) should handle association if implemented correctly.
    # Wait, setBFInterval! in OctreeBuilder.jl was not fully shown/checked. 
    # Let's assume it works for now or check it later.
    
    # 5. Translation
    # Run translation for all levels
    for iLevel in 2:octree.nLevels
        translate!(octree.levels[iLevel])
    end
    
    # Check disaggG
    @test !isempty(level2.disaggG)
    
    # 6. Disaggregation
    # Downward pass
    for iLevel in 2:(octree.nLevels - 1)
        parentLevel = octree.levels[iLevel]
        childLevel = octree.levels[iLevel + 1]
        disaggregate_downward!(parentLevel, childLevel)
    end
    
    # Leaf disaggregation
    ZI = zeros(ComplexF64, num_basis(basis))
    disaggregate_leaf!(leafLevel, basis, operator, ZI, sorted_ids)
    
    # Check if ZI has been updated (might be small or zero if no incident field/translation was effective)
    # But at least it should run without error.
    @test length(ZI) == num_basis(basis)

end

@testset "MLFMA Octree Construction (Random Points)" begin
    # Generate random points
    n_points = 1000
    points = rand(3, n_points) .* 10.0 # 10x10x10 box
    
    leafCubeEdgel = 1.0
    
    octree, sorted_ids = build_octree(points, leafCubeEdgel)
    
    @test octree isa OctreeInfo
    @test octree.nLevels > 0
    @test length(sorted_ids) == n_points
    
    # Check levels
    for (id, level) in octree.levels
        @test level isa LevelInfo
        @test level.ID == id
        @test length(level.cubes) > 0
    end
    
    # Check leaf level
    leafLevel = octree.levels[octree.nLevels]
    @test leafLevel.isleaf
    
    # Check top level (level 2)
    if octree.nLevels >= 2
        topLevel = octree.levels[2]
        @test !topLevel.isleaf
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 14.6: MLFMAOperatorMPI
# ─────────────────────────────────────────────────────────────────────────────
@testset "MLFMAOperatorMPI (single-process)" begin
    # Initialize MPI for single-process use
    if !MPI.Initialized()
        MPI.Init()
    end
    comm = MPI.COMM_WORLD

    # Small sphere mesh for fast tests
    radius = 0.5
    freq   = 300e6   # λ ≈ 1 m
    mesh   = generate_sphere_mesh(radius, 6, 12)
    basis  = RWGBasis(mesh)
    N      = num_basis(basis)
    @test N > 0

    efie          = EFIE(freq)
    leafCubeEdgel = 0.3

    # Serial reference
    mlfma_serial = MLFMAOperator(efie, basis, leafCubeEdgel)
    @test mlfma_serial isa MLFMAOperator

    # Distributed (single process — should be identical to serial)
    mlfma_mpi = MLFMAOperatorMPI(efie, basis, leafCubeEdgel; comm = comm)
    @test mlfma_mpi isa MLFMAOperatorMPI
    @test size(mlfma_mpi) == size(mlfma_serial)

    # Random test vector
    import Random: seed!
    seed!(1234)
    x = randn(ComplexF64, N)

    y_serial = similar(x)
    y_mpi    = similar(x)
    mul!(y_serial, mlfma_serial, x)
    mul!(y_mpi,    mlfma_mpi,    x)

    # On P=1, results must be bitwise identical
    @test maximum(abs.(y_mpi .- y_serial)) == 0.0
    @test norm(y_mpi) ≈ norm(y_serial)
end

@testset "MLFMAOperator preserves physical basis ordering" begin
    radius = 0.5
    freq = 300e6
    mesh = generate_sphere_mesh(radius, 6, 12)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)

    set_frequency!(freq)
    Z_direct = assemble_impedance_matrix(efie, basis)
    source = PlaneWave(freq, π / 2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)
    I_direct = Z_direct \ V

    op = MLFMAOperator(efie, basis, 0.3)
    precond = LUPreconditioner(lu(op.Z_near))

    solver_phys = GMRESSolver(restart = 20, maxiter = 80, tol = 1e-3, verbose = false)
    I_phys = solve!(solver_phys, op, V; Pl = precond)

    solver_sorted = GMRESSolver(restart = 20, maxiter = 80, tol = 1e-3, verbose = false)
    I_sorted_rhs = solve!(solver_sorted, op, V[op.sorted_ids]; Pl = precond)

    err_phys = norm(I_phys - I_direct) / norm(I_direct)
    err_sorted_rhs = norm(I_sorted_rhs - I_direct) / norm(I_direct)

    @test err_phys < 1e-8
    @test err_sorted_rhs > 1e-1
end

@testset "issue #20: interpolationCSCMatCal honours clamped per-direction order" begin
    Interp = EMMoMSuite.FastAlgorithms.MLFMA.Interpolation

    # A parent level can have fewer directional samples than the default interpolation
    # order (nInterp = 6) when the octree is deep (fine leaf boxes). In that case the code
    # clamps the order per-direction (nlocalInterpTheta / nlocalInterpPhi), but it used to
    # keep sizing/indexing the sparse arrays with the *un-clamped* nlocalInterp, which
    # crashed with a BoundsError in the θ step (or a sparse length-mismatch in φ).
    parentL, parentPoles = Interp.levelIntegralInfoCal(0.02; λ = 1.0)  # L=4 -> nθ=5 (<6), nφ=10
    childL,  childPoles  = Interp.levelIntegralInfoCal(0.01; λ = 1.0)  # L=3 -> nθ=4, nφ=8

    @test parentL ≤ 4                                   # θ order is clamped (< 6)
    @test length(parentPoles.Xθs) < 6
    info = Interp.interpolationCSCMatCal(parentPoles, childPoles, 6)   # must not throw
    @test info isa Interp.LagrangeInterpInfo

    # φ-clamped path (issue #20's reported symptom): too-fine parent with nφ < 6.
    phiParentL, phiParentPoles = Interp.levelIntegralInfoCal(0.001; λ = 1.0)
    phiChildL,  phiChildPoles  = Interp.levelIntegralInfoCal(0.0005; λ = 1.0)
    @test length(phiParentPoles.Xϕs) < 6
    info2 = Interp.interpolationCSCMatCal(phiParentPoles, phiChildPoles, 6)
    @test info2 isa Interp.LagrangeInterpInfo
end

@testset "issue #22: MLFMAOperator element cache + scratch reuse" begin
    freq = 300e6
    mesh = generate_sphere_mesh(0.5, 6, 12)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)

    op = MLFMAOperator(efie, basis, 0.3)

    # Element geometry/quadrature must be precomputed once at construction
    @test op.element_cache !== nothing
    infos, gqs = op.element_cache
    @test length(infos) == 1 && length(gqs) == 1
    @test infos[1] !== nothing && gqs[1] !== nothing

    # Scratch-buffer reuse must not corrupt results: repeated matvecs are
    # deterministic, and mul! into a caller buffer agrees with `*`.
    x = ones(ComplexF64, size(op, 1))
    y1 = op * x
    y2 = op * x
    @test y1 == y2
    y3 = similar(x)
    mul!(y3, op, x)
    @test y3 == y1
    @test all(isfinite, y1)
end

@testset "issue #22: BlockJacobiPreconditioner on MLFMAOperator" begin
    freq = 300e6
    mesh = generate_sphere_mesh(0.5, 6, 12)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    op = MLFMAOperator(efie, basis, 0.3)

    P = BlockJacobiPreconditioner(op)
    x = ones(ComplexF64, size(op, 1))
    y = P \ x
    @test length(y) == length(x)
    @test all(isfinite, y)
end

@testset "issue #22 #3: adaptive tree-uniform near_range" begin
    OB = EMMoMSuite.FastAlgorithms.MLFMA.OctreeBuilder

    # Calibration: GD2V gate fixture (λ=1, w=0.1 → L=9) must reproduce
    # the empirically passing near_range=7 (kR_min = 8·kw ≈ 5.03 ≥ 0.55·L = 4.95)
    @test OB.adaptive_near_range(1.0, 0.1, 9) == 7
    # issue-#22 bench leaf (λ=0.3, w=0.06 → L=12): kR_min = 6·kw ≈ 7.54 ≥ 6.6
    @test OB.adaptive_near_range(0.3, 0.06, 12) == 5
    # floor at 1 for electrically large cubes
    @test OB.adaptive_near_range(1.0, 2.0, 9) == 1

    freq = 300e6
    mesh = generate_sphere_mesh(0.5, 6, 12)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    centers = reduce(hcat, [bf.center for bf in basis.functions])

    # Default (adaptive, leaf-derived tree-uniform radius): every level's
    # smallest far pair must satisfy kR_min ≥ 0.55·L — correct M2L by
    # construction (problem 3). nr = 7 > 2 also triggers the efficiency
    # guidance warning.
    octree, _ = @test_logs (:warn, r"Adaptive near_range = 7") match_mode = :any begin
        build_octree(centers, 0.1; λ = 1.0)
    end
    k = 2π
    leaf_nr = OB.adaptive_near_range(
        1.0,
        octree.levels[octree.nLevels].cubeEdgel,
        octree.levels[octree.nLevels].L,
    )
    for levelID = 2:octree.nLevels
        level = octree.levels[levelID]
        kR_min = Inf
        n_far = 0
        for cube in level.cubes, fid in cube.farneighbors
            far = level.cubes[fid]
            off = (cube.ID3D[1] - far.ID3D[1],
                   cube.ID3D[2] - far.ID3D[2],
                   cube.ID3D[3] - far.ID3D[3])
            kR_min = min(kR_min, k * norm(off) * level.cubeEdgel)
            n_far += 1
        end
        @test kR_min ≥ 0.55 * level.L - 1e-9
    end

    # Near/far TILING: no existing cube pair may fall through both lists.
    # For every ordered leaf pair (i, j): offset d ≤ nr ⇒ j ∈ neighbors(i);
    # nr < d ≤ 2nr+1 ⇒ j ∈ farneighbors(i). (d > 2nr+1 is handled at coarser
    # levels — that is the hierarchical tiling.) This is the regression that
    # catches per-level radii breaking the parent-window invariant.
    leaf = octree.levels[octree.nLevels]
    for i in 1:leaf.nCubes
        ci = leaf.cubes[i]
        for j in 1:leaf.nCubes
            j == i && continue
            cj = leaf.cubes[j]
            d = max(abs(ci.ID3D[1] - cj.ID3D[1]),
                    abs(ci.ID3D[2] - cj.ID3D[2]),
                    abs(ci.ID3D[3] - cj.ID3D[3]))
            d ≤ leaf_nr || d ≤ 2 * leaf_nr + 1 || continue
            if d ≤ leaf_nr
                @test (j in ci.neighbors)
            else
                @test (j in ci.farneighbors)
            end
        end
    end

    # Explicit near_range is still honored: leaf neighbors within ±near_range
    oct2, _ = build_octree(centers, 0.1; λ = 1.0, near_range = 1)
    leaf = oct2.levels[oct2.nLevels]
    maxoff = 0
    for cube in leaf.cubes, n in cube.neighbors
        off = max(abs(cube.ID3D[1] - leaf.cubes[n].ID3D[1]),
                  abs(cube.ID3D[2] - leaf.cubes[n].ID3D[2]),
                  abs(cube.ID3D[3] - leaf.cubes[n].ID3D[3]))
        maxoff = max(maxoff, off)
    end
    @test maxoff ≤ 1
end

@testset "issue #22 #3: scaled M2L stabilization — documented limitation" begin
    # The scaled diagonalization (Ergül & Karaosmanoğlu, URSI GA 2014)
    # stabilizes the h_l series for SUBWAVELENGTH boxes (kw ≪ 1 — the paper's
    # λ/250 demo). At electrical-size boxes (kw ≈ 2π·0.3 here) the short-range
    # error is dominated by evanescent content that NEITHER the unscaled NOR
    # the s^l-scaled real-spectrum M2L captures — empirically the scaled form
    # is orders of magnitude WORSE (rel ≈ 1e4 vs direct). The correct remedy
    # at conventional frequencies is the adaptive near_range (GD2V-calibrated);
    # the multipole-domain (spectral) reformulation was investigated and
    # FALSIFIED at integration (see the SphericalHarmonics testset below and
    # docs/dev/m2l_short_range_spectral.md §6).
    # This testset pins the documented, opt-in behavior: the operator runs and
    # produces finite output; it makes no accuracy claim.
    freq = 300e6
    mesh = generate_sphere_mesh(0.5, 6, 12)
    basis = RWGBasis(mesh)
    efie = EFIE(freq)
    set_frequency!(freq)
    x = ones(ComplexF64, num_basis(basis))

    Z = assemble_impedance_matrix(efie, basis)
    y_dense = Z * x

    op_s = MLFMAOperator(efie, basis, 0.3, Val(:Lagrange2Step), 1;
                         m2l_stabilization = :scaled)
    y_s = op_s * x
    rel_s = norm(y_s - y_dense) / norm(y_dense)
    @info "scaled M2L at nr=1 (documented limitation: not accuracy-recovering " *
          "at electrical-size boxes)" rel_s

    @test all(isfinite, y_s)
    @test all(isfinite, y_dense)
end

# ─────────────────────────────────────────────────────────────────────────────
# issue #22 问题 3 调研产物：SphericalHarmonics 数学工具（单元自检）
#
# 背景（负结果，2025 阶段 3 集成结论）：谱域（多极子域）M2L——分析
# `F̂ = YᵀW·F` → 稠密 Gaunt 平移矩阵 `M` → 综合 `G = const·W·Y·M·F̂`——
# 在理想化原型（带限 F、共享求积）中成立（scripts/verification/
# verify_m2l_spectral_shortrange.jl 自检全绿），但在真实管线中被证伪：
# 对角方案的逐点采样求积在 kR < L（h_l 放大区，实际叶层常态）依赖 GL
# 采样对截断级数尾部的广义求和/混叠相消才收敛于真值（GD2X nr=2 对角
# rel=3.6e-3）；带限往返破坏该相消（谱域 rel ≥ 3.6e2，"精确 Gaunt 配对"
# 的低阶谱系数比真值大 ~10³×）。机制与数据：docs/dev/
# m2l_short_range_spectral.md §6。本模块保留为经机器精度验证的数学
# 工具（球谐分析/综合、Gaunt 表、谱域平移矩阵——供分析与后续研究）。
# ─────────────────────────────────────────────────────────────────────────────

@testset "issue #22 #3: SphericalHarmonics 单元自检" begin
    SH = EMMoMSuite.FastAlgorithms.MLFMA.SphericalHarmonics
    Interp = EMMoMSuite.FastAlgorithms.MLFMA.Interpolation

    # 1) GL 极点网格上球谐矩阵正交归一（degree ≤ L 精确 ⇒ 机器精度）
    for L in (9, 12)
        Lc, poles = Interp.levelIntegralInfoCal(0.1; λ = 0.5, L_min = L)
        Y, YtW = SH.sh_matrix(poles, L)
        S = YtW * Y
        @test maximum(abs, S - I) < 1e-12
    end

    # 2) 超额带宽公式：与解析公式一致、单调、有下限
    @test SH.spectral_bandwidth(1.2566; digits = 2.0) == 5   # GD2V 叶盒 kw
    @test SH.spectral_bandwidth(2.5133; digits = 2.0) == 9   # 0.2λ 盒
    @test SH.spectral_bandwidth(0.0) == 0

    # 3) M 矩阵单谐解析锚点：M_{a,(0,0)}·√(4π) = 4π(-j)^{l_a} h_{l_a} conj(Y_a(R̂))
    L = 9
    k = 4π
    Rvec = [0.2, 0.35, 0.91]
    R = norm(Rvec); R̂ = Rvec / R; kR = k * R
    θR = acos(clamp(R̂[3], -1, 1)); φR = atan(R̂[2], R̂[1])
    Mb = SH.m2l_matrix(L, Rvec, k)
    plgndr = SH.plgndr
    normf(l, m) = sqrt((2l + 1) / (4π)) *
        exp(0.5 * (SpecialFunctions.loggamma(l - abs(m) + 1) -
                   SpecialFunctions.loggamma(l + abs(m) + 1)))
    err = 0.0
    a00 = 1   # idx(0,0) = 0²+(0+0+1) = 1
    for l = 0:L, m = -l:l
        ai = l^2 + (m + l + 1)
        θpart = m >= 0 ? normf(l, m) * plgndr(l, m, cos(θR)) :
                (-1.0)^(-m) * normf(l, m) * plgndr(l, -m, cos(θR))
        Ya = θpart * exp(im * m * φR)
        rhs = 4π * (-im)^l *
              (SpecialFunctions.sphericalbesselj(l, kR) -
               im * SpecialFunctions.sphericalbessely(l, kR)) * conj(Ya)
        err = max(err, abs(Mb[ai, a00] * sqrt(4π) - rhs) / max(abs(rhs), 1e-30))
    end
    @test err < 1e-10

    # 4) 行列对称截断：Leff 截断只保留 l ≤ Leff 的行×列子块
    Mb6 = SH.m2l_matrix(L, Rvec, k; Leff = 6)
    @test maximum(abs, Mb6[1:49, 1:49] .- Mb[1:49, 1:49]) == 0  # (l ≤ 6)² 子块不变
    @test maximum(abs, Mb6[:, 50:end]) == 0                     # l_b > 6 列全零
    @test maximum(abs, Mb6[50:end, :]) == 0                     # l_a > 6 行全零
end

