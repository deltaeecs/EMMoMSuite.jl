using Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))
using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.FastAlgorithms.MLFMA
using EMSuite.FastAlgorithms.MLFMA.MLFMAOperatorModule: compute_interaction_term
using LinearAlgebra
using StaticArrays
using Printf
using Random

function debug_scfie_performance()
    println("==================================================")
    println("   Debug: SCFIE Interaction Performance           ")
    println("==================================================")

    # 1. Parameters
    freq = 1.2e9
    EMSuite.Utilities.Parameters.set_frequency!(freq)
    
    # 2. Mesh
    mesh_file = joinpath(@__DIR__, "../../../MoM_AllinOne/meshfiles/plate_and_metal_1dot2GHz.nas")
    if !isfile(mesh_file)
        println("Error: Mesh file not found: $mesh_file")
        return
    end
    
    println("Loading mesh...")
    surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file)
    
    # 3. Basis
    println("Setting up Bases...")
    basis_surf = RWGBasis(surf_mesh)
    basis_vol = SWGBasis(vol_mesh)
    
    # 4. Operator
    eps_r = 2.0 * (1.0 - 0.0002im)
    permittivities = fill(eps_r, num_elements(vol_mesh))
    scfie = SCFIE(freq, permittivities; alpha=0.0)
    
    # Pre-create sub-operators
    efie_op = EFIE(freq)
    mfie_op = MFIE(freq)
    vefie_op = VEFIE(freq, permittivities)
    
    # 5. Select Random Elements for Interaction
    println("Benchmarking Interaction Terms...")
    
    n_trials = 10000
    
    # Surface-Surface
    println("\n--- Surface-Surface (EFIE/MFIE) ---")
    tri_test = get_triangle_info(surf_mesh, basis_surf, 1)
    tri_src = get_triangle_info(surf_mesh, basis_surf, 2)
    
    # Warmup
    EMSuite.FastAlgorithms.MLFMA.MLFMAOperatorModule.compute_interaction_rwg_rwg(scfie, tri_test, tri_src, 1, 1, efie_op, mfie_op)
    
    t_ss = @elapsed for i in 1:n_trials
        EMSuite.FastAlgorithms.MLFMA.MLFMAOperatorModule.compute_interaction_rwg_rwg(scfie, tri_test, tri_src, 1, 1, efie_op, mfie_op)
    end
    println("Time for $n_trials interactions: $(t_ss) s")
    println("Average time per interaction: $(t_ss/n_trials*1e6) μs")
    
    # Volume-Volume
    println("\n--- Volume-Volume (VEFIE) ---")
    tet_infos = get_tetrahedra_info(vol_mesh, basis_vol, permittivities)
    tet_test = tet_infos[1]
    tet_src = tet_infos[2]
    
    # Warmup
    EMSuite.FastAlgorithms.MLFMA.MLFMAOperatorModule.compute_interaction_swg_swg(scfie, tet_test, tet_src, 1, 1, vefie_op)
    
    t_vv = @elapsed for i in 1:n_trials
        EMSuite.FastAlgorithms.MLFMA.MLFMAOperatorModule.compute_interaction_swg_swg(scfie, tet_test, tet_src, 1, 1, vefie_op)
    end
    println("Time for $n_trials interactions: $(t_vv) s")
    println("Average time per interaction: $(t_vv/n_trials*1e6) μs")
    
    # Surface-Volume
    println("\n--- Surface-Volume (Coupling) ---")
    
    # Warmup
    EMSuite.FastAlgorithms.MLFMA.MLFMAOperatorModule.compute_interaction_rwg_swg(scfie, tri_test, tet_src, 1, 1)
    
    t_sv = @elapsed for i in 1:n_trials
        EMSuite.FastAlgorithms.MLFMA.MLFMAOperatorModule.compute_interaction_rwg_swg(scfie, tri_test, tet_src, 1, 1)
    end
    println("Time for $n_trials interactions: $(t_sv) s")
    println("Average time per interaction: $(t_sv/n_trials*1e6) μs")
    
    # Volume-Volume (Cached)
    println("\n--- Volume-Volume (VEFIE) Cached ---")
    
    # Warmup
    Z_ts, _ = EMSuite.IntegralEquations.VEFIEModule.vefie_element_interaction(vefie_op, tet_test, tet_src)
    val = Z_ts[1, 1]
    
    t_vv_cached = @elapsed for i in 1:n_trials
        # Simulate loop over 16 basis pairs
        # Compute matrix ONCE
        Z_ts, _ = EMSuite.IntegralEquations.VEFIEModule.vefie_element_interaction(vefie_op, tet_test, tet_src)
        # Access 16 times
        for m in 1:4, n in 1:4
            val = Z_ts[m, n]
        end
    end
    # t_vv_cached is for n_trials * 16 interactions
    avg_time_cached = t_vv_cached / (n_trials * 16) * 1e6
    println("Time for $(n_trials*16) interactions: $(t_vv_cached) s")
    println("Average time per interaction (Cached): $(avg_time_cached) μs")
    
    # Surface-Volume (Cached)
    println("\n--- Surface-Volume (Coupling) Cached ---")
    t_sv_cached = @elapsed for i in 1:n_trials
        # Simulate loop over 12 basis pairs (3x4)
        Z_sv, _ = EMSuite.IntegralEquations.SCFIEModule.scfie_coupling_interaction(scfie, tri_test, tet_src)
        for m in 1:3, n in 1:4
            val = Z_sv[m, n]
        end
    end
    avg_time_sv_cached = t_sv_cached / (n_trials * 12) * 1e6
    println("Time for $(n_trials*12) interactions: $(t_sv_cached) s")
    println("Average time per interaction (Cached): $(avg_time_sv_cached) μs")
    
    # Estimate Total Time
    nnz_est = 24.6e6
    # Weighted average? Assume mostly Volume-Volume
    avg_time_est = avg_time_cached
    total_est = avg_time_est * 1e-6 * nnz_est
    println("\n--- Estimation (Cached) ---")
    println("Estimated Total Assembly Time (1 Thread): $(total_est) s")
    
end

debug_scfie_performance()
