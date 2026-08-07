"""
    完整精度效率基准测试
    运行方式: julia --project=EMMoMSuite -t4 benchmark/run_full_benchmark.jl
    从 MoM/ 根目录运行
"""
using EMMoMSuite
using LinearAlgebra
using Printf
using Statistics

const BASE_DIR = joinpath(@__DIR__, "..")      # EMMoMSuite/
const MOM_DIR = joinpath(BASE_DIR, "..")       # MoM/
const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne")
const MOM_BASICS_DIR = joinpath(@__DIR__, "..", "deps", "fixtures", "Basics")
const LEGACY_BASELINE_DIR = joinpath(BASE_DIR, "test_results", "legacy_baseline")

# LU preconditioner wrapper (needed for IterativeSolvers ldiv! interface)
struct LUPreconditioner
    F
end
LinearAlgebra.ldiv!(y, P::LUPreconditioner, x) = (y .= P.F \ x)
LinearAlgebra.ldiv!(P::LUPreconditioner, x) = (x .= P.F \ x)

# Try loading CSV for comparison
has_csv = try
    @eval using CSV, DataFrames
    true
catch
    false
end

println("=" ^ 80)
println("EMMoMSuite 精度效率基准测试")
println("Julia $(VERSION), Threads: $(Threads.nthreads())")
println("=" ^ 80)

# ============================================================
# Test 1: SEFIE Direct (Jet, 100MHz, EFIE, LU)
# ============================================================
function test_sefie_direct()
    println("\n" * "=" ^ 60)
    println("Test 1: SEFIE Direct (Jet, 100MHz)")
    println("=" ^ 60)

    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
    if !isfile(mesh_file)
        println("  SKIP: mesh file not found: $mesh_file")
        return nothing
    end

    mesh = read_nas_mesh(mesh_file, scale=1.0)
    freq = 1e8
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("  Unknowns: $N, Elements: $(num_elements(mesh))")

    efie = EFIE(freq)

    # Assembly
    t_asm = @elapsed Z = assemble_impedance_matrix(efie, basis)
    println("  Assembly: $(round(t_asm, digits=2))s")

    # Excitation (Legacy convention: direction -x, polarization +z)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)

    # Solve
    t_slv = @elapsed I_coeff = solve!(LUSolver(), Z, V)
    println("  Solve: $(round(t_slv, digits=2))s")
    println("  Total: $(round(t_asm + t_slv, digits=2))s")

    # RCS
    θs = collect(LinRange(-π, π, 721))
    ϕs = [0.0, π/2]
    RCS_res = radarCrossSection(θs, ϕs, I_coeff, basis)
    RCS_linear = RCS_res[2]  # Total RCS, linear (m²)
    RCS_dBsm = 10 * log10.(RCS_linear)

    # Compare with Legacy
    baseline_file = joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv")
    if has_csv && isfile(baseline_file)
        df = CSV.read(baseline_file, DataFrame)
        d0 = RCS_dBsm[:, 1] .- df.RCS_Phi0_dBsm
        d90 = RCS_dBsm[:, 2] .- df.RCS_Phi90_dBsm

        println("  --- Phi=0 vs Legacy ---")
        println("  Mean Diff: $(round(mean(d0), digits=3)) dB")
        println("  RMSE: $(round(sqrt(mean(d0.^2)), digits=3)) dB")
        println("  Max |Diff|: $(round(maximum(abs.(d0)), digits=3)) dB")

        println("  --- Phi=90 vs Legacy ---")
        println("  Mean Diff: $(round(mean(d90), digits=3)) dB")
        println("  RMSE: $(round(sqrt(mean(d90.^2)), digits=3)) dB")
        println("  Max |Diff|: $(round(maximum(abs.(d90)), digits=3)) dB")

        println("  --- 典型角度 (Phi=0) ---")
        for (i, th) in enumerate(θs)
            td = round(th * 180 / π, digits=1)
            if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
                @printf("  θ=%+7.1f° | EMMoMSuite=%+8.2f | Legacy=%+8.2f | Diff=%+6.2f dB\n",
                    td, RCS_dBsm[i,1], df.RCS_Phi0_dBsm[i], d0[i])
            end
        end

        return (test="SEFIE Direct", N=N, t_asm=t_asm, t_slv=t_slv, t_total=t_asm+t_slv,
                mean_d0=mean(d0), rmse_d0=sqrt(mean(d0.^2)),
                mean_d90=mean(d90), rmse_d90=sqrt(mean(d90.^2)))
    else
        println("  WARNING: Legacy baseline not found (CSV)")
        return (test="SEFIE Direct", N=N, t_asm=t_asm, t_slv=t_slv, t_total=t_asm+t_slv)
    end
end

# ============================================================
# Test 2: SEFIE MLFMA (Jet, 100MHz, EFIE, GMRES)
# ============================================================
function test_sefie_mlfma()
    println("\n" * "=" ^ 60)
    println("Test 2: SEFIE MLFMA (Jet, 100MHz)")
    println("=" ^ 60)

    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
    if !isfile(mesh_file)
        println("  SKIP: mesh file not found")
        return nothing
    end

    mesh = read_nas_mesh(mesh_file, scale=1.0)
    freq = 1e8
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    λ = 299792458.0 / freq
    leaf_size = 0.35 * λ
    println("  Unknowns: $N, λ=$(round(λ, digits=3))m")

    efie = EFIE(freq)

    # MLFMAOperator handles octree, near-field internally. Sorting is transparent.
    t_setup = @elapsed mlfma_op = MLFMAOperator(efie, basis, leaf_size)
    println("  MLFMA Setup: $(round(t_setup, digits=2))s")
    println("  Octree Levels: $(mlfma_op.octree.nLevels)")

    # Excitation (no sorting needed — MLFMA is transparent)
    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(efie, source, basis)

    # Preconditioner
    P = LUPreconditioner(lu(mlfma_op.Z_near))

    # Solve
    solver = GMRESSolver(restart=50, maxiter=100, tol=1e-3, verbose=true)
    t_slv = @elapsed I_coeff = solve!(solver, mlfma_op, V, Pl=P)
    println("  GMRES Solve: $(round(t_slv, digits=2))s")
    println("  Total: $(round(t_setup + t_slv, digits=2))s")

    # RCS (no unsort needed)
    θs = collect(LinRange(-π, π, 721))
    ϕs = [0.0, π/2]
    RCS_res = radarCrossSection(θs, ϕs, I_coeff, basis)
    RCS_linear = RCS_res[2]
    RCS_dBsm = 10 * log10.(RCS_linear)

    # Compare with Legacy
    baseline_file = joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv")
    if has_csv && isfile(baseline_file)
        df = CSV.read(baseline_file, DataFrame)
        d0 = RCS_dBsm[:, 1] .- df.RCS_Phi0_dBsm

        println("  --- vs Legacy Direct (Phi=0) ---")
        println("  Mean Diff: $(round(mean(d0), digits=3)) dB")
        println("  RMSE: $(round(sqrt(mean(d0.^2)), digits=3)) dB")
        println("  Max |Diff|: $(round(maximum(abs.(d0)), digits=3)) dB")

        println("  --- 典型角度 ---")
        for (i, th) in enumerate(θs)
            td = round(th * 180 / π, digits=1)
            if td in [-180.0, -90.0, 0.0, 90.0, 180.0]
                @printf("  θ=%+7.1f° | MLFMA=%+8.2f | Legacy=%+8.2f | Diff=%+6.2f dB\n",
                    td, RCS_dBsm[i,1], df.RCS_Phi0_dBsm[i], d0[i])
            end
        end

        return (test="SEFIE MLFMA", N=N, t_setup=t_setup, t_slv=t_slv, t_total=t_setup+t_slv,
                mean_d0=mean(d0), rmse_d0=sqrt(mean(d0.^2)))
    end
    return (test="SEFIE MLFMA", N=N, t_setup=t_setup, t_slv=t_slv, t_total=t_setup+t_slv)
end

# ============================================================
# Test 3: SCFIE Direct (Jet, 100MHz, CFIE α=0.5, LU)
# ============================================================
function test_scfie_direct()
    println("\n" * "=" ^ 60)
    println("Test 3: SCFIE Direct (Jet, 100MHz, CFIE α=0.5)")
    println("=" ^ 60)

    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
    if !isfile(mesh_file)
        println("  SKIP: mesh file not found")
        return nothing
    end

    mesh = read_nas_mesh(mesh_file, scale=1.0)
    freq = 1e8
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    println("  Unknowns: $N")

    cfie = CFIE(freq, 0.5)

    t_asm = @elapsed Z = assemble_impedance_matrix(cfie, basis)
    println("  Assembly: $(round(t_asm, digits=2))s")

    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)

    t_slv = @elapsed I_coeff = solve!(LUSolver(), Z, V)
    println("  Solve: $(round(t_slv, digits=2))s")
    println("  Total: $(round(t_asm + t_slv, digits=2))s")

    # RCS
    θs = collect(LinRange(-π, π, 721))
    ϕs = [0.0, π/2]
    RCS_res = radarCrossSection(θs, ϕs, I_coeff, basis)
    RCS_linear = RCS_res[2]
    RCS_dBsm = 10 * log10.(RCS_linear)

    # Compare with Legacy EFIE Direct
    baseline_file = joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv")
    if has_csv && isfile(baseline_file)
        df = CSV.read(baseline_file, DataFrame)
        d0 = RCS_dBsm[:, 1] .- df.RCS_Phi0_dBsm

        println("  --- vs Legacy EFIE Direct (Phi=0) ---")
        println("  Mean Diff: $(round(mean(d0), digits=3)) dB")
        println("  RMSE: $(round(sqrt(mean(d0.^2)), digits=3)) dB")

        println("  --- 典型角度 ---")
        for (i, th) in enumerate(θs)
            td = round(th * 180 / π, digits=1)
            if td in [-180.0, 0.0, 180.0]
                @printf("  θ=%+7.1f° | CFIE=%+8.2f | Legacy EFIE=%+8.2f | Diff=%+6.2f dB\n",
                    td, RCS_dBsm[i,1], df.RCS_Phi0_dBsm[i], d0[i])
            end
        end

        return (test="SCFIE Direct", N=N, t_asm=t_asm, t_slv=t_slv, t_total=t_asm+t_slv,
                mean_d0=mean(d0), rmse_d0=sqrt(mean(d0.^2)))
    end
    return (test="SCFIE Direct", N=N, t_asm=t_asm, t_slv=t_slv, t_total=t_asm+t_slv)
end

# ============================================================
# Test 4: SCFIE MLFMA MatVec 精度 (TriTetra, mixed freq)
# ============================================================
function test_scfie_mlfma_matvec()
    println("\n" * "=" ^ 60)
    println("Test 4: SCFIE MLFMA MatVec (TriTetra)")
    println("=" ^ 60)

    mesh_file = joinpath(MOM_BASICS_DIR, "meshfiles", "TriTetra.nas")
    if !isfile(mesh_file)
        println("  SKIP: mesh file not found: $mesh_file")
        return nothing
    end

    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file, scale=0.001)

    eps_r_val = 2.0 * (1 - 0.001im)

    results = Dict()

    for (label, freq) in [("2GHz (pure near)", 2e9), ("4GHz (mixed)", 4e9)]
        set_frequency!(freq)
        λ = 299792458.0 / freq
        perms = fill(eps_r_val, num_elements(vol_mesh))

        surf_basis = RWGBasis(surf_mesh)
        vol_basis = SWGBasis(vol_mesh)
        N_s = num_basis(surf_basis)
        N_v = num_basis(vol_basis)
        N = N_s + N_v

        scfie = SCFIE(freq, perms; alpha=0.5)

        # Direct matrix
        Z_direct = assemble_impedance_matrix(scfie, surf_basis, vol_basis)

        # MLFMA
        leaf_size = 0.25 * λ
        mlfma_op = MLFMAOperator(scfie, [surf_basis, vol_basis], leaf_size)

        # MatVec comparison (no sorting needed — transparent)
        x = randn(ComplexF64, N)
        y_direct = Z_direct * x
        y_mlfma = mlfma_op * x

        rel_err = norm(y_direct - y_mlfma) / norm(y_direct)

        println("  [$label] N=$N (S=$N_s, V=$N_v), Rel Err: $(round(rel_err*100, digits=4))%")
        results[label] = rel_err
    end

    return (test="SCFIE MLFMA MatVec", results=results)
end

# ============================================================
# Test 5: EFIE MLFMA MatVec (TriTetra surface, 4GHz)
# ============================================================
function test_efie_mlfma_matvec()
    println("\n" * "=" ^ 60)
    println("Test 5: EFIE MLFMA MatVec (TriTetra surf, 4GHz)")
    println("=" ^ 60)

    mesh_file = joinpath(MOM_BASICS_DIR, "meshfiles", "TriTetra.nas")
    if !isfile(mesh_file)
        println("  SKIP: mesh not found")
        return nothing
    end

    surf_mesh, _ = read_mixed_nas_mesh(mesh_file, scale=0.001)
    freq = 4e9
    set_frequency!(freq)
    λ = 299792458.0 / freq
    basis = RWGBasis(surf_mesh)
    N = num_basis(basis)
    println("  Unknowns: $N")

    efie = EFIE(freq)

    Z_direct = assemble_impedance_matrix(efie, basis)

    leaf_size = 0.25 * λ
    mlfma_op = MLFMAOperator(efie, basis, leaf_size)

    x = randn(ComplexF64, N)
    y_direct = Z_direct * x
    y_mlfma = mlfma_op * x

    rel_err = norm(y_direct - y_mlfma) / norm(y_direct)
    println("  MatVec Relative Error: $(round(rel_err*100, digits=4))%")

    return (test="EFIE MLFMA MatVec", N=N, err=rel_err)
end

# ============================================================
# Test 6: VEFIE MLFMA MatVec (TriTetra vol, mixed freq)
# ============================================================
function test_vefie_mlfma_matvec()
    println("\n" * "=" ^ 60)
    println("Test 6: VEFIE MLFMA MatVec (TriTetra vol)")
    println("=" ^ 60)

    mesh_file = joinpath(MOM_BASICS_DIR, "meshfiles", "TriTetra.nas")
    if !isfile(mesh_file)
        println("  SKIP: mesh not found")
        return nothing
    end

    # Also test pure tetra mesh
    mesh_tetra = joinpath(MOM_BASICS_DIR, "meshfiles", "Tetra.nas")
    
    results = Dict()

    # TriTetra volume part at 4GHz
    _, vol_mesh = read_mixed_nas_mesh(mesh_file, scale=0.001)
    freq = 4e9
    set_frequency!(freq)
    λ = 299792458.0 / freq
    basis = SWGBasis(vol_mesh)
    N = num_basis(basis)
    eps_r = fill(complex(2.0, -0.002), num_elements(vol_mesh))
    println("  [TriTetra 4GHz] Unknowns: $N")

    vefie = VEFIE(freq, eps_r)
    Z_direct = assemble_impedance_matrix(vefie, basis)

    leaf_size = 0.25 * λ
    mlfma_op = MLFMAOperator(vefie, basis, leaf_size)

    x = randn(ComplexF64, N)
    y_direct = Z_direct * x
    y_mlfma = mlfma_op * x
    rel_err = norm(y_direct - y_mlfma) / norm(y_direct)
    println("  MatVec Relative Error: $(round(rel_err*100, digits=4))%")
    results["TriTetra 4GHz"] = rel_err

    # Pure Tetra at 1.5GHz
    if isfile(mesh_tetra)
        tet_mesh = read_nas_mesh(mesh_tetra, scale=0.001)
        freq2 = 1.5e9
        set_frequency!(freq2)
        λ2 = 299792458.0 / freq2
        basis2 = SWGBasis(tet_mesh)
        N2 = num_basis(basis2)
        eps_r2 = fill(complex(2.0, -0.002), num_elements(tet_mesh))
        println("  [Tetra 1.5GHz] Unknowns: $N2")

        vefie2 = VEFIE(freq2, eps_r2)
        Z_d2 = assemble_impedance_matrix(vefie2, basis2)
        leaf2 = 0.25 * λ2
        op2 = MLFMAOperator(vefie2, basis2, leaf2)

        x2 = randn(ComplexF64, N2)
        y_d2 = Z_d2 * x2
        y_m2 = op2 * x2
        err2 = norm(y_d2 - y_m2) / norm(y_d2)
        println("  MatVec Relative Error: $(round(err2*100, digits=4))%")
        results["Tetra 1.5GHz"] = err2
    end

    return (test="VEFIE MLFMA MatVec", results=results)
end

# ============================================================
# Test 7: SCFIE MLFMA Full Solve (Jet, 100MHz, CFIE α=0.5)
# ============================================================
function test_scfie_mlfma()
    println("\n" * "=" ^ 60)
    println("Test 7: SCFIE MLFMA (Jet, 100MHz, CFIE α=0.5)")
    println("=" ^ 60)

    mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
    if !isfile(mesh_file)
        println("  SKIP: mesh file not found")
        return nothing
    end

    mesh = read_nas_mesh(mesh_file, scale=1.0)
    freq = 1e8
    set_frequency!(freq)
    basis = RWGBasis(mesh)
    N = num_basis(basis)
    λ = 299792458.0 / freq
    leaf_size = 0.35 * λ
    println("  Unknowns: $N")

    cfie = CFIE(freq, 0.5)

    t_setup = @elapsed mlfma_op = MLFMAOperator(cfie, basis, leaf_size)
    println("  MLFMA Setup: $(round(t_setup, digits=2))s")
    println("  Octree Levels: $(mlfma_op.octree.nLevels)")

    source = PlaneWave(freq, π/2, π, [0.0, 0.0, 1.0])
    V = excitation_vector(cfie, source, basis)

    P = LUPreconditioner(lu(mlfma_op.Z_near))
    solver = GMRESSolver(restart=50, maxiter=200, tol=1e-3, verbose=true)
    t_slv = @elapsed I_coeff = solve!(solver, mlfma_op, V, Pl=P)
    println("  GMRES Solve: $(round(t_slv, digits=2))s")
    println("  Total: $(round(t_setup + t_slv, digits=2))s")

    # RCS
    θs = collect(LinRange(-π, π, 721))
    ϕs = [0.0, π/2]
    RCS_res = radarCrossSection(θs, ϕs, I_coeff, basis)
    RCS_linear = RCS_res[2]
    RCS_dBsm = 10 * log10.(RCS_linear)

    # Compare with Legacy EFIE Direct
    baseline_file = joinpath(LEGACY_BASELINE_DIR, "SEFIE_Direct_Jet.csv")
    if has_csv && isfile(baseline_file)
        df = CSV.read(baseline_file, DataFrame)
        d0 = RCS_dBsm[:, 1] .- df.RCS_Phi0_dBsm

        println("  --- vs Legacy EFIE Direct (Phi=0) ---")
        println("  Mean Diff: $(round(mean(d0), digits=3)) dB")
        println("  RMSE: $(round(sqrt(mean(d0.^2)), digits=3)) dB")

        return (test="SCFIE MLFMA", N=N, t_setup=t_setup, t_slv=t_slv, t_total=t_setup+t_slv,
                mean_d0=mean(d0), rmse_d0=sqrt(mean(d0.^2)))
    end
    return (test="SCFIE MLFMA", N=N, t_setup=t_setup, t_slv=t_slv, t_total=t_setup+t_slv)
end

# ============================================================
# Test 8: SCFIE 耦合互易性 (TriTetra, 300MHz)
# ============================================================
function test_scfie_reciprocity()
    println("\n" * "=" ^ 60)
    println("Test 8: SCFIE 耦合互易性 (TriTetra, 300MHz)")
    println("=" ^ 60)

    mesh_file = joinpath(MOM_BASICS_DIR, "meshfiles", "TriTetra.nas")
    if !isfile(mesh_file)
        println("  SKIP")
        return nothing
    end

    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file, scale=0.001)
    freq = 3e8
    set_frequency!(freq)
    eps_r = fill(2.0 + 0im, num_elements(vol_mesh))

    surf_basis = RWGBasis(surf_mesh)
    vol_basis = SWGBasis(vol_mesh)
    N_s = num_basis(surf_basis)
    N_v = num_basis(vol_basis)
    println("  RWG: $N_s, SWG: $N_v")

    scfie = SCFIE(freq, eps_r; alpha=0.5)
    Z = assemble_impedance_matrix(scfie, surf_basis, vol_basis)

    Z_SS = Z[1:N_s, 1:N_s]
    Z_SV = Z[1:N_s, N_s+1:end]
    Z_VS = Z[N_s+1:end, 1:N_s]
    Z_VV = Z[N_s+1:end, N_s+1:end]

    println("  ||Z_SS||: $(round(norm(Z_SS), sigdigits=4))")
    println("  ||Z_SV||: $(round(norm(Z_SV), sigdigits=4))")
    println("  ||Z_VS||: $(round(norm(Z_VS), sigdigits=4))")
    println("  ||Z_VV||: $(round(norm(Z_VV), sigdigits=4))")

    return (test="SCFIE Reciprocity", N_s=N_s, N_v=N_v)
end

# ============================================================
# Run All Tests
# ============================================================
all_results = []

for (label, fn) in [
    ("Test 1", test_sefie_direct),
    ("Test 2", test_sefie_mlfma),
    ("Test 3", test_scfie_direct),
    ("Test 4", test_scfie_mlfma_matvec),
    ("Test 5", test_efie_mlfma_matvec),
    ("Test 6", test_vefie_mlfma_matvec),
    ("Test 7", test_scfie_mlfma),
    ("Test 8", test_scfie_reciprocity),
]
    try
        r = fn()
        push!(all_results, r)
    catch ex
        println("  ERROR in $label: $ex")
        for line in stacktrace(catch_backtrace())[1:min(10, end)]
            println("    ", line)
        end
        push!(all_results, nothing)
    end
end

# ============================================================
# Summary
# ============================================================
println("\n" * "=" ^ 80)
println("BENCHMARK SUMMARY")
println("=" ^ 80)

for r in all_results
    r === nothing && continue
    print("  $(r.test) | N=$(hasproperty(r, :N) ? r.N : "N/A")")
    if hasproperty(r, :t_total)
        @printf(" | Time=%.1fs", r.t_total)
    end
    if hasproperty(r, :rmse_d0)
        @printf(" | RMSE=%.3f dB", r.rmse_d0)
    end
    if hasproperty(r, :err)
        @printf(" | MatVec Err=%.4f%%", r.err*100)
    end
    println()
end

println("\nDONE")
