# run_case_mpi.jl — MPI multi-process case runner
#
# Runs each case with a DISTRIBUTED MLFMA + distributed GMRES solve
# (DistributedSPAIPreconditioner, BlockJacobi fallback) so the shipped report
# reflects an MPI run. Surface cases use RWG; volume (VEFIE) cases use SWG
# with a per-object leaf strategy leaf = clamp(Lmax/6, 0.2λ, 0.5λ).
# A dense-LU solution on rank 0 serves as the cross-method reference.
#
# Usage:  mpiexec -n <P> julia --project=benchmark benchmark/cases/run_case_mpi.jl [CASE...]
# e.g.:   mpiexec -n 2 julia --project=benchmark benchmark/cases/run_case_mpi.jl G1_sphere_CFIE_Mie
#
# Rank 0 writes the usual artifact set (report.md, rcs.csv, perf.csv, PNGs,
# spec.toml, case.log); every rank participates in the distributed solve.

using MPI
using EMMoMSuite
using LinearAlgebra, Statistics, Printf
import LinearAlgebra: ldiv!
using EMMoMSuite.FastAlgorithms.MLFMA: MLFMAOperatorMPI
using EMMoMSuite.Parallel:
    DistributedSPAIPreconditioner, DistributedBlockJacobiPreconditioner,
    apply_mpi_preconditioner!

include(joinpath(@__DIR__, "CaseRunner.jl"))
using .CaseRunner
using .CaseRunner: postprocess_rcs, postprocess_farfield, reference_rcs,
                   build_operator, check_orientation!, CaseResult,
                   _write_rcs_csv, _plot_rcs, _plot_farfield, _plot_current_views,
                   _plot_geometry_and_mesh, _publication_report, _deg,
                   build_volume_operator,
                   _write_spec_snapshot, RESULT_ROOT, GEO_DIR

# apply_mpi_preconditioner! 适配成 IterativeSolvers 左预条件 ldiv! 接口
struct FuncPrecond{F}
    f::F
end
ldiv!(y, P::FuncPrecond, x) = P.f(y, x)
function ldiv!(P::FuncPrecond, x)
    y = similar(x)
    P.f(y, x)
    copyto!(x, y)
    return x
end

# 分布式 VEFIE 体积分支 — 镜像 CaseRunner._run_volume_case：
# gmsh 体网格 → 区域材料绑定 → 体算子（各 rank 相同）→ SWG 基 → 激励
# → 分布式 MLFMA + GMRES → rank 0 稠密 LU 交叉验证 + 报告。
#
# leaf 策略：leaf = clamp(物体最长边/6, 0.2λ, 0.5λ)，避免过小 leaf 触发
# near_range>2 警告与树退化。主 leaf 用于正式产物；另取 clamp 边界的另一端
# 做一次对比求解（GMRES 残差 + RCS RMSE vs LU），环境变量 EMMOM_LEAF_COMPARE=0
# 可跳过；EMMOM_MPI_LEAF=<米> 可覆盖主 leaf。
function _run_volume_case_mpi(spec, comm, rank, P, outdir)
    # ---- 1. geometry -> gmsh tetrahedral mesh + Physical Volume regions ---
    geo_file = isabspath(spec.geo) ? spec.geo : joinpath(GEO_DIR, spec.geo)
    t0 = time()
    mesh, region_tags = generate_gmsh_volume(geo_file; mesh_size = spec.mesh_size)
    t_mesh = time() - t0

    # ---- 1b. bind materials by region name and validate -------------------
    bm = bind_regions(mesh, region_tags,
                      Dict(r.surface => r.material for r in spec.regions))
    validate_bindings(bm) ||
        error("case $(spec.name): unbound tetrahedra (Physical Volume missing?)")
    tetnum = mesh.tetnum
    nnodes = size(mesh.node, 2)
    rank == 0 && @printf("  mesh: %d tetrahedra, %d nodes, regions = %s (%.1f s)\n",
                         tetnum, nnodes, join(keys(region_tags), ", "), t_mesh)

    set_frequency!(spec.freq)
    source = PlaneWave(spec.freq, spec.theta_inc, spec.phi_inc, spec.pol)

    # ---- 2. volume operator (identical on every rank) + SWG basis ---------
    op, solve_mesh, ctx = build_volume_operator(bm, spec.regions, spec.freq;
                                                ie = spec.ie)
    ctx.boundary_fallback &&
        error("case $(spec.name): MPI volume path requires dielectric SWG " *
              "(all-PEC regions fall back to RWG — use the surface path)")
    swg = SWGBasis(solve_mesh)
    nbasis = num_basis(swg)
    V = excitation_vector(op, source, swg, ctx.permittivities)  # perms: 每四面体 εᵣ
    rank == 0 && println("  assembly ($(nameof(typeof(op))), $nbasis unknowns)")

    # ---- 2b. leaf strategy ------------------------------------------------
    λ = 3e8 / spec.freq
    ext = vec(maximum(mesh.node; dims = 2) .- minimum(mesh.node; dims = 2))
    Lmax = maximum(ext)                       # 物体最长边 [m]
    leaf_main = clamp(Lmax / 6, 0.2λ, 0.5λ)
    haskey(ENV, "EMMOM_MPI_LEAF") && (leaf_main = parse(Float64, ENV["EMMOM_MPI_LEAF"]))
    leaf_alt = leaf_main ≈ 0.5λ ? 0.2λ : 0.5λ
    rank == 0 && @printf("  object longest edge %.3f m (λ = %.3f m): leaf = %.4f m (= %.2fλ), comparison leaf = %.4f m\n",
                         Lmax, λ, leaf_main, leaf_main / λ, leaf_alt)

    # ---- 3. distributed MLFMA + preconditioned GMRES ----------------------
    t0 = time()
    op_mpi = MLFMAOperatorMPI(op, swg, leaf_main; comm = comm)
    t_assembly = time() - t0
    rank == 0 && @printf("  distributed MLFMA near-field assembly (leaf %.4f m): %.1f s\n",
                         leaf_main, t_assembly)

    gmres_tol = parse(Float64, get(ENV, "EMMOM_GMRES_TOL", "1e-3"))
    pname, I, t_solve, rel_res = _mpi_gmres_volume(op_mpi, V, comm, rank, P,
                                                   gmres_tol, leaf_main)

    # ---- 3b. leaf comparison (one extra solve, rank timing only) ----------
    if get(ENV, "EMMOM_LEAF_COMPARE", "1") != "0"
        t0 = time()
        op_mpi2 = MLFMAOperatorMPI(op, swg, leaf_alt; comm = comm)
        _, _, t_solve2, rel_res2 = _mpi_gmres_volume(op_mpi2, V, comm, rank, P,
                                                     gmres_tol, leaf_alt)
        t2 = time() - t0
        rank == 0 && @printf("  leaf compare: leaf %.4f m → rel-res %.2e in %.1f s | leaf %.4f m → rel-res %.2e in %.1f s\n",
                             leaf_main, rel_res, t_solve, leaf_alt, rel_res2, t2)
        op_mpi2 = nothing; GC.gc()
    end

    if rank == 0
        # ---- 4. RCS -------------------------------------------------------
        θa = collect(Float64.(range(0.0, pi; length = spec.n_theta)))
        ϕs = collect(spec.phi_cuts)
        t0 = time()
        rcs_dB = postprocess_rcs(ctx, θa, ϕs, I, swg)      # [nθ, nϕ] dBsm
        t_rcs = time() - t0

        # ---- 4b. cross-check: dense LU on rank 0 --------------------------
        mlfma_ok = false
        mlfma_rmse = Float64[]
        t_lu = NaN
        rcs_lu = nothing
        try
            t0 = time()
            Z = assemble_impedance_matrix(op, swg)
            I_lu = Z \ V
            t_lu = time() - t0
            rcs_lu = postprocess_rcs(ctx, θa, ϕs, I_lu, swg)
            mlfma_rmse = [sqrt(mean((rcs_dB[:, j] .- rcs_lu[:, j]).^2))
                          for j in eachindex(ϕs)]
            mlfma_ok = true
            for (j, φ) in enumerate(ϕs)
                @printf("  MPI GMRES vs dense LU, phi=%6.1f°: RMSE = %.3f dB (LU %.1f s)\n",
                        _deg(φ), mlfma_rmse[j], t_lu)
            end
        catch err
            @warn "dense-LU cross-check failed" spec.name err
        end

        # ---- 5. far field + plots + report --------------------------------
        t0 = time()
        FF = postprocess_farfield(ctx, θa, ϕs, I, swg, source, nbasis)
        # rcs_dB 是分布式 MLFMA 主解；把稠密 LU 作为对比曲线/列写入
        nan2 = fill(NaN, length(θa), length(ϕs))
        lu_dB = mlfma_ok ? rcs_lu : nothing
        _write_rcs_csv(spec, outdir, θa, ϕs, rcs_dB, nan2, false; mlfma_dB = lu_dB)
        _plot_rcs(spec, outdir, θa, ϕs, rcs_dB, nan2, false;
                  mlfma_dB = lu_dB, main_label = "MLFMA (MPI)",
                  mlfma_label = "MoM (dense LU)")
        _plot_farfield(spec, outdir, θa, ϕs, FF)
        t_plots = time() - t0
        @printf("  plots: %.1f s\n", t_plots)

        res = CaseResult(spec, outdir, mesh, tetnum, nnodes, nbasis,
                         t_mesh, t_assembly, t_solve, t_rcs, t_plots,
                         false, Float64[], mlfma_ok, mlfma_rmse, t_lu)
        _publication_report(res, rcs_dB, fill(NaN, length(θa), length(ϕs)), false;
                            geometry = :volume, region_tags = region_tags, ctx = ctx,
                            mpi = "MPI distributed GMRES (P = $P ranks, " *
                                  "$pname precond, leaf = $(round(leaf_main; digits = 4)) m, " *
                                  "rel-res = $(round(rel_res; sigdigits = 2)))")
        println("  report: ", joinpath(outdir, "report.md"))
    end
    return nothing
end

# SPAI 优先；MethodError 或收敛差（rel-res > 10×tol）时回退 BlockJacobi
function _mpi_gmres_volume(op_mpi, V, comm, rank, P, tol, leaf)
    pname = "SAI"
    P_mpi = try
        DistributedSPAIPreconditioner(op_mpi)
    catch err
        err isa MethodError || rethrow()
        rank == 0 && @warn "DistributedSPAIPreconditioner not available for this " *
                           "operator; falling back to DistributedBlockJacobiPreconditioner" err
        nothing
    end
    I = nothing
    rel_res = Inf
    for attempt in 1:2
        if P_mpi === nothing
            pname = "BlockJacobi"
            P_mpi = DistributedBlockJacobiPreconditioner(op_mpi)
        end
        pwrap = FuncPrecond((y, t) -> apply_mpi_preconditioner!(y, P_mpi, t))
        solver = GMRESSolver(restart = 200, maxiter = 400, tol = tol,
                             verbose = rank == 0)
        t0 = time()
        I = solve!(solver, op_mpi, V; Pl = pwrap)
        t_solve = time() - t0
        # 分布式相对残差 ‖Ax−b‖/‖b‖（matvec 后各 rank 均持有全量 y）
        r = mul!(similar(V), op_mpi, I) .- V
        g = MPI.Allreduce([real(dot(r, r)), real(dot(V, V))], +, comm)
        rel_res = sqrt(g[1] / g[2])
        rank == 0 && @printf("  distributed GMRES (%s, P = %d, leaf %.4f m, tol = %.0e): %.1f s, rel-res = %.2e\n",
                             pname, P, leaf, tol, t_solve, rel_res)
        (rel_res <= 10 * tol || pname == "BlockJacobi") && return pname, I, t_solve, rel_res
        rank == 0 && @warn "poor GMRES convergence with SPAI; retrying with BlockJacobi" rel_res
        P_mpi = nothing
    end
    return pname, I, t_solve, rel_res
end

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    names = isempty(ARGS) ? String[] : ARGS
    specs = CaseRunner.load_cases()
    isempty(names) || (specs = filter(s -> s.name in names, specs))
    # MPI 体积分支（VEFIE/SWG）；表面积分走原有 RWG 路径
    isvol(spec) = !isempty(spec.regions)

    for spec in specs
        rank == 0 && println("=== Case $(spec.name) (MPI, P=$P) ===")
        outdir = joinpath(RESULT_ROOT, spec.name)
        if rank == 0
            mkpath(outdir)
            _write_spec_snapshot(spec, outdir)
        end

        if isvol(spec)
            _run_volume_case_mpi(spec, comm, rank, P, outdir)
            MPI.Barrier(comm)
            continue
        end

        # ---- 1. geometry -> gmsh mesh (每个 rank 独立生成同一网格) --------
        geo_file = isabspath(spec.geo) ? spec.geo : joinpath(GEO_DIR, spec.geo)
        t0 = time()
        mesh = generate_gmsh_from_file(geo_file; mesh_size = spec.mesh_size, dim = spec.dim)
        t_mesh = time() - t0
        check_orientation!(spec, mesh)
        trinum = mesh.trinum
        nnodes = size(mesh.node, 2)
        basis = RWGBasis(mesh)
        nbasis = num_basis(basis)
        rank == 0 && @printf("  mesh: %d triangles, %d nodes, %d RWG unknowns (%.1f s)\n",
                             trinum, nnodes, nbasis, t_mesh)

        set_frequency!(spec.freq)
        source = PlaneWave(spec.freq, spec.theta_inc, spec.phi_inc, spec.pol)
        op, ctx = build_operator(spec.interfaces, spec.freq; ie = spec.ie,
                                 alpha = spec.alpha)
        V = excitation_vector(op, source, basis)

        # ---- 2-3. distributed MLFMA assembly + GMRES solve ----------------
        λ = 3e8 / spec.freq
        t0 = time()
        op_mpi = op isa PMCHW ?
                 PMCHWMLFMAOperatorMPI(op, basis, λ / 2; comm = comm) :
                 MLFMAOperatorMPI(op, basis, λ / 2; comm = comm)
        t_assembly = time() - t0
        rank == 0 && @printf("  distributed MLFMA near-field assembly: %.1f s\n", t_assembly)

        t0 = time()
        if op isa PMCHW
            P_mpi = DistributedBlockJacobiPreconditioner(op_mpi)
        else
            P_mpi = DistributedSPAIPreconditioner(op_mpi)
        end
        pwrap = FuncPrecond((y, t) -> apply_mpi_preconditioner!(y, P_mpi, t))
        # GMRES 收敛阈值：默认 1e-3；可用环境变量 EMMOM_GMRES_TOL 覆盖
        gmres_tol = parse(Float64, get(ENV, "EMMOM_GMRES_TOL", "1e-3"))
        solver = GMRESSolver(restart = 200, maxiter = 400, tol = gmres_tol,
                             verbose = rank == 0)
        I = solve!(solver, op_mpi, V; Pl = pwrap)
        t_solve = time() - t0
        rank == 0 && @printf("  distributed GMRES solve (SAI, P = %d, tol = %.0e): %.1f s\n",
                             P, gmres_tol, t_solve)
        MPI.Barrier(comm)

        if rank == 0
            # ---- 4. RCS + Mie reference ----------------------------------
            t0 = time()
            θa = collect(Float64.(range(0.0, pi; length = spec.n_theta)))
            ϕs = collect(spec.phi_cuts)
            rcs_dB = postprocess_rcs(ctx, θa, ϕs, I, basis)   # [nθ, nϕ] dBsm

            mie_dB = Matrix{Float64}(undef, length(θa), length(ϕs)); fill!(mie_dB, NaN)
            mie_ok = false
            if spec.mie_radius !== nothing
                try
                    for (j, φ) in enumerate(ϕs)
                        ref = reference_rcs(ctx, spec.mie_radius, spec.freq,
                                            θa, φ, spec.theta_inc, spec.phi_inc, spec.pol)
                        mie_dB[:, j] .= ref[2]
                    end
                    mie_ok = true
                catch err
                    @warn "Mie reference failed" spec.name err
                end
            end
            t_rcs = time() - t0

            rmse = [mie_ok ? sqrt(mean((rcs_dB[:, j] .- mie_dB[:, j]).^2)) : NaN
                    for j in eachindex(ϕs)]
            for (j, φ) in enumerate(ϕs)
                mie_ok && @printf("  RCS vs Mie, phi=%6.1f°: RMSE = %.3f dB\n",
                                  _deg(φ), rmse[j])
            end

            # ---- 4b. cross-check: dense LU on rank 0 ----------------------
            mlfma_ok = false
            mlfma_rmse = Float64[]
            t_lu = NaN
            rcs_lu = nothing
            try
                t0 = time()
                Z = assemble_impedance_matrix(op, basis)
                I_lu = Z \ V
                t_lu = time() - t0
                rcs_lu = postprocess_rcs(ctx, θa, ϕs, I_lu, basis)
                mlfma_rmse = [sqrt(mean((rcs_dB[:, j] .- rcs_lu[:, j]).^2))
                              for j in eachindex(ϕs)]
                mlfma_ok = true
                for (j, φ) in enumerate(ϕs)
                    @printf("  MPI GMRES vs dense LU, phi=%6.1f°: RMSE = %.3f dB (LU %.1f s)\n",
                            _deg(φ), mlfma_rmse[j], t_lu)
                end
            catch err
                @warn "dense-LU cross-check failed" spec.name err
            end

            # ---- 5. far field + plots -------------------------------------
            t0 = time()
            FF = postprocess_farfield(ctx, θa, ϕs, I, basis, source, nbasis)
            # rcs_dB 是分布式 MLFMA 主解；把稠密 LU 作为对比曲线/列写入
            lu_dB = mlfma_ok ? rcs_lu : nothing
            _write_rcs_csv(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok;
                           mlfma_dB = lu_dB)
            _plot_rcs(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok;
                      mlfma_dB = lu_dB, main_label = "MLFMA (MPI)",
                      mlfma_label = "MoM (dense LU)")
            _plot_farfield(spec, outdir, θa, ϕs, FF)
            try   # surface-current distribution (J part only; PMCHW M ignored here)
                I_j = ctx.layout === :jm ? I[1:nbasis] : I
                J = geoElectricJCal(I_j, basis)
                _plot_current_views(spec, mesh, J, outdir)
            catch err
                @warn "surface-current plot failed" spec.name err
            end
            _plot_geometry_and_mesh(spec, mesh, outdir)
            t_plots = time() - t0
            @printf("  plots: %.1f s\n", t_plots)

            res = CaseResult(spec, outdir, mesh, trinum, nnodes, nbasis,
                             t_mesh, t_assembly, t_solve, t_rcs, t_plots,
                             mie_ok, rmse, mlfma_ok, mlfma_rmse, t_lu)
            _publication_report(res, rcs_dB, mie_dB, mie_ok;
                                geometry = :surface,
                                mpi = "MPI distributed GMRES (P = $P ranks)")
            println("  report: ", joinpath(outdir, "report.md"))
        end
        MPI.Barrier(comm)
    end
    MPI.Finalize()
    return nothing
end

main()
