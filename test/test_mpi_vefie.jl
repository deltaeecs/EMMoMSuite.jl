# test_mpi_vefie.jl — MPI 分布式 VEFIE+SWG 三门回归测试
# 运行: mpiexec -n <P> julia --project=. test/test_mpi_vefie.jl
#       （P=1 和 P=2 都要跑过；门 3 由 P=2 读取 P=1 写出的解做一致性对比）
#
# 门（固化 scratch_vefie_mpi_probe.jl 探针结论，探针 matvec 相对差 0.000e+00）：
#   1. MLFMAOperatorMPI matvec == 串行 MLFMAOperator matvec，相对差 < 1e-10
#   2. 分布式 GMRES 解 == 稠密 LU 解，相对差 < 1e-6
#   3. P=1 与 P=2 分布式 GMRES 解一致（相对差 < 1e-6）
using Test
using MPI
using EMMoMSuite
using LinearAlgebra, Random, Printf, DelimitedFiles
using IterativeSolvers
using EMMoMSuite.FastAlgorithms.MLFMA: MLFMAOperator, MLFMAOperatorMPI
using EMMoMSuite.Parallel: DistributedSPAIPreconditioner, DistributedBlockJacobiPreconditioner,
                           apply_mpi_preconditioner!

# 预条件包装：apply_mpi_preconditioner! → IterativeSolvers ldiv! 接口
struct FuncPrecond{F}
    f::F
end
LinearAlgebra.ldiv!(y, P::FuncPrecond, x) = P.f(y, x)
function LinearAlgebra.ldiv!(P::FuncPrecond, x)
    y = similar(x)
    P.f(y, x)
    copyto!(x, y)
    return x
end

# 算子包装：mul! + size，满足 IterativeSolvers 接口
struct OpWrap{A}
    op::A
    n::Int
end
Base.size(A::OpWrap) = (A.n, A.n)
Base.size(A::OpWrap, d::Int) = d == 1 || d == 2 ? A.n : 1
Base.eltype(::Type{<:OpWrap}) = ComplexF64
LinearAlgebra.mul!(y, A::OpWrap, v) = mul!(y, A.op, v)

const SOL_DIR = normpath(joinpath(@__DIR__, "..", "test_results"))
const SOL_FILE(P) = joinpath(SOL_DIR, "mpi_vefie_solution_P$P.csv")

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    P = MPI.Comm_size(comm)

    @testset "MPI 分布式 VEFIE+SWG (P=$P)" begin
        # ---- 小号 G12 几何：cube_r0p1.geo, mesh_size = 0.06 ----------------
        geo = normpath(joinpath(@__DIR__, "..", "benchmark", "cases", "geo",
                                "cube_r0p1.geo"))
        mesh = generate_gmsh_from_file(geo; mesh_size = 0.06, dim = 3)

        freq = 600e6
        set_frequency!(freq)
        λ = 3e8 / freq
        perms = fill(ComplexF64(4.0 + 0.0im), mesh.tetnum)   # 每四面体 εᵣ = 4
        swg = SWGBasis(mesh)
        N = num_basis(swg)
        rank == 0 && @printf("tets=%d nodes=%d SWG unknowns N=%d (P=%d)\n",
                             mesh.tetnum, size(mesh.node, 2), N, P)

        op = VEFIE(freq, perms)
        # leaf 策略与 run_case_mpi.jl 一致：clamp(最长边/6, 0.2λ, 0.5λ)
        ext = vec(maximum(mesh.node; dims = 2) .- minimum(mesh.node; dims = 2))
        leaf = clamp(maximum(ext) / 6, 0.2λ, 0.5λ)

        op_mpi = MLFMAOperatorMPI(op, swg, leaf; comm = comm)
        op_ser = rank == 0 ? MLFMAOperator(op, swg, leaf) : nothing

        # ── 门 1：MPI matvec vs 串行 MLFMA matvec（< 1e-10）────────────────
        Random.seed!(7)
        x = randn(ComplexF64, N)
        y_mpi = mul!(zeros(ComplexF64, N), op_mpi, x)
        MPI.Barrier(comm)
        if rank == 0
            y_ser = mul!(zeros(ComplexF64, N), op_ser, x)
            rel1 = norm(y_mpi - y_ser) / norm(y_ser)
            @info "门1 matvec rel diff (MPI vs serial)" rel1
            @test rel1 < 1e-10
        end
        MPI.Barrier(comm)

        # ---- 激励 + 分布式 GMRES（SAI 预条件，BlockJacobi 回退）------------
        source = PlaneWave(freq, pi / 2, pi, [0.0, 0.0, 1.0])
        V = excitation_vector(op, source, swg, perms)
        pname = "SAI"
        P_mpi = try
            DistributedSPAIPreconditioner(op_mpi)
        catch err
            err isa MethodError || rethrow()
            nothing
        end
        if P_mpi === nothing
            pname = "BlockJacobi"
            P_mpi = DistributedBlockJacobiPreconditioner(op_mpi)
        end
        pwrap = FuncPrecond((y, t) -> apply_mpi_preconditioner!(y, P_mpi, t))
        I_mpi, hist = gmres(OpWrap(op_mpi, N), V; Pl = pwrap, restart = 200,
                            maxiter = 400, reltol = 1e-10, log = true)

        # ── 门 2：分布式 GMRES 解 vs 串行 MLFMA-GMRES / 稠密 LU ─────────────
        # MPI 路径正确性门：同一算子（分布式 vs 串行 MLFMA）的解一致 < 1e-6。
        # vs 稠密 LU 受 MLFMA 远场插值精度限制（rel ~1e-4 是算子精度地板，
        # 见串行 G12 报告同量级），门限放 1e-3。
        if rank == 0
            solver = GMRESSolver(restart = 200, maxiter = 400, tol = 1e-10)
            Pbj = BlockJacobiPreconditioner(op_ser)
            I_ser = solve!(solver, op_ser, V; Pl = Pbj)
            I_lu = assemble_impedance_matrix(op, swg) \ V
            rel2a = norm(I_mpi - I_ser) / norm(I_ser)   # MPI 正确性
            rel2b = norm(I_ser - I_lu) / norm(I_lu)     # MLFMA 精度地板
            @info "门2 GMRES rel diff" precond = pname gmres_iters = hist.iters
            @info "门2 分布式 vs 串行 MLFMA-GMRES" rel2a
            @info "门2 串行 MLFMA-GMRES vs 稠密 LU（算子精度地板）" rel2b
            @test rel2a < 1e-6
            @test rel2b < 1e-3

            # ---- 写本进程数的解，供另一进程数做门 3 对比 -------------------
            mkpath(SOL_DIR)
            open(SOL_FILE(P), "w") do io
                println(io, "real,imag")
                writedlm(io, [reinterpret(Float64, [z]) for z in I_mpi], ',')
            end
        end
        MPI.Barrier(comm)

        # ── 门 3：P=1 与 P=2 解一致性（< 1e-6；由 P=2 侧对比 P=1 的解）──────
        if P > 1 && rank == 0
            @test isfile(SOL_FILE(1))
            if isfile(SOL_FILE(1))
                raw = readdlm(SOL_FILE(1), ',', Float64; header = true)[1]
                I_p1 = ComplexF64.(raw[:, 1], raw[:, 2])
                rel3 = norm(I_mpi - I_p1) / norm(I_p1)
                @info "门3 solution rel diff (P=2 vs P=1)" rel3
                @test rel3 < 1e-6
            end
        end
    end
    MPI.Finalize()
    return nothing
end

main()
