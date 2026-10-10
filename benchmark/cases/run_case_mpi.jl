# run_case_mpi.jl — MPI multi-process case runner
#
# Runs a surface case with a DISTRIBUTED MLFMA + distributed GMRES solve
# (DistributedSPAIPreconditioner) so the shipped report reflects an MPI run.
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

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)
    names = isempty(ARGS) ? String[] : ARGS
    specs = CaseRunner.load_cases()
    isempty(names) || (specs = filter(s -> s.name in names, specs))
    # MPI 路径只支持表面积分（体积 VEFIE 尚无分布式算子，另行处理）
    specs = filter(s -> isempty(s.regions), specs)

    for spec in specs
        rank == 0 && println("=== Case $(spec.name) (MPI, P=$P) ===")
        outdir = joinpath(RESULT_ROOT, spec.name)
        if rank == 0
            mkpath(outdir)
            _write_spec_snapshot(spec, outdir)
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
        solver = GMRESSolver(restart = 200, maxiter = 400, tol = 1e-6,
                             verbose = rank == 0)
        I = solve!(solver, op_mpi, V; Pl = pwrap)
        t_solve = time() - t0
        rank == 0 && @printf("  distributed GMRES solve (SAI, P = %d): %.1f s\n", P, t_solve)
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
            _write_rcs_csv(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok)
            _plot_rcs(spec, outdir, θa, ϕs, rcs_dB, mie_dB, mie_ok)
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
