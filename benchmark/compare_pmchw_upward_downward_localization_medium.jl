using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using LinearAlgebra
using Random
using Printf
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule: PMCHWMLFMAOperator, PMCHWMLFMAErrorBudget, aggregate_leaf_pmchw!
using EMSuite.FastAlgorithms.MLFMA.PMCHWMLFMAOperatorModule
using EMSuite.FastAlgorithms.MLFMA.Aggregation: aggregate_upward!
using EMSuite.FastAlgorithms.MLFMA.Disaggregation: disaggregate_downward!
using EMSuite.FastAlgorithms.MLFMA.Translation: translate!
using EMSuite.IntegralEquations.PMCHWModule: assemble_K_pmchw_offdiag
using EMSuite.IntegralEquations.Impedance: get_triangles_info
using EMSuite.Geometry: GaussQuadratureInfo

function make_fixture()
    mesh = generate_sphere_mesh(0.5, 10, 20)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(300e6, 4.0)
    return pmchw, basis
end

function make_probe(N, which::Symbol; seed = 43)
    Random.seed!(seed)
    x = zeros(ComplexF64, 2N)
    if which === :M
        x[(N + 1):(2N)] .= randn(ComplexF64, N)
    elseif which === :J
        x[1:N] .= randn(ComplexF64, N)
    else
        error("unsupported probe selector: $which")
    end
    x ./= norm(x)
    return x
end

function clear_agg!(oct)
    for (_, lv) in oct.levels
        isdefined(lv, :aggS) && fill!(lv.aggS, zero(eltype(lv.aggS)))
        isdefined(lv, :disaggG) && fill!(lv.disaggG, zero(eltype(lv.disaggG)))
    end
end

function near_pairs_basis(octree, sorted_ids)
    leaf_level = octree.levels[octree.nLevels]
    pairs = Set{Tuple{Int,Int}}()

    for i_cube in 1:leaf_level.nCubes
        cube = leaf_level.cubes[i_cube]
        isempty(cube.bfInterval) && continue
        test_ids = [sorted_ids[s] for s in cube.bfInterval]

        for neigh_idx in cube.neighbors
            neigh_cube = leaf_level.cubes[neigh_idx]
            isempty(neigh_cube.bfInterval) && continue
            src_ids = [sorted_ids[s] for s in neigh_cube.bfInterval]
            for i in test_ids, j in src_ids
                push!(pairs, (i, j))
            end
        end
    end

    return pairs
end

function direct_m_pass_refs(pmchw, basis, op)
    k1_real = abs(real(pmchw.k1))
    eta1_real = abs(real(pmchw.eta1))
    em_k0 = assemble_K_pmchw_offdiag(basis, pmchw.k0)
    em_k1 = assemble_K_pmchw_offdiag(basis, k1_real)
    hm_k0 = assemble_impedance_matrix(efie_from_keta(pmchw.k0, pmchw.eta0, im * pmchw.k0 / (pmchw.eta0 * 16pi)), basis)
    hm_k1 = assemble_impedance_matrix(efie_from_keta(k1_real, eta1_real, im * pmchw.k1 / (pmchw.eta1 * 16pi)), basis)

    near_k0 = near_pairs_basis(op.octree0, op.sorted_ids0)
    near_k1 = near_pairs_basis(op.octree1, op.sorted_ids1)
    for (i, j) in near_k0
        em_k0[i, j] = 0
        hm_k0[i, j] = 0
    end
    for (i, j) in near_k1
        em_k1[i, j] = 0
        hm_k1[i, j] = 0
    end

    return (
        k0 = (EM = em_k0, HM = hm_k0),
        k1 = (EM = em_k1, HM = hm_k1),
    )
end

function exact_level_agg(level, basis, x, sorted_ids, x_range, k; quad_order = 4)
    offset = first(x_range) - 1
    FT = eltype(level.cubeEdgel)
    CT = Complex{FT}
    JK = CT(im * k)

    poles = level.poles
    nPoles = length(poles.r̂sθsϕs)
    exact = zeros(CT, nPoles, 2, level.nCubes)

    tri_info = get_triangles_info(basis.mesh, basis)
    gq = GaussQuadratureInfo(:Triangle, quad_order, FT)
    n_qp = length(gq.weight)

    poles_r̂ = [p.r̂ for p in poles.r̂sθsϕs]
    poles_θhat = [p.θhat for p in poles.r̂sθsϕs]
    poles_ϕhat = [p.ϕhat for p in poles.r̂sθsϕs]

    for iCube in 1:level.nCubes
        cube = level.cubes[iCube]
        cubeCenter = cube.center

        for bfID_sorted in cube.bfInterval
            bfID_orig = sorted_ids[bfID_sorted]
            coeff = x[bfID_orig + offset]
            abs(coeff) < 1e-12 && continue

            bf = basis.functions[bfID_orig]
            for i_supp in 1:2
                tri_idx = bf.support[i_supp]
                tri_idx == 0 && continue

                tri = tri_info[tri_idx]
                v_all = tri.vertices
                local_edge = bf.local_edge_idx[i_supp]
                v_opp = v_all[:, local_edge]
                sign_supp = bf.signs[i_supp]

                for i_qp in 1:n_qp
                    L = gq.coordinate[:, i_qp]
                    r = v_all * L
                    rho = r - v_opp
                    r_local = r - cubeCenter
                    factor_vec = sign_supp * bf.edge_length / 2 * gq.weight[i_qp] * coeff

                    for iPole in 1:nPoles
                        phase = exp(JK * dot(poles_r̂[iPole], r_local))
                        vec = rho * factor_vec * phase
                        exact[iPole, 1, iCube] += dot(poles_θhat[iPole], vec)
                        exact[iPole, 2, iCube] += dot(poles_ϕhat[iPole], vec)
                    end
                end
            end
        end
    end

    return exact
end

function worst_cube_metrics(approx, exact)
    worst_cube = 0
    worst_rel = -1.0
    worst_abs = 0.0
    worst_norm = 0.0

    for iCube in axes(exact, 3)
        approx_cube = @view approx[:, :, iCube]
        exact_cube = @view exact[:, :, iCube]
        exact_norm = norm(exact_cube)
        abs_err = norm(approx_cube - exact_cube)
        rel_err = abs_err / (exact_norm + 1e-30)
        if rel_err > worst_rel
            worst_rel = rel_err
            worst_abs = abs_err
            worst_norm = exact_norm
            worst_cube = iCube
        end
    end

    return (cube = worst_cube, rel = worst_rel, abs = worst_abs, ref = worst_norm)
end

function capture_disagg_snapshots(octree)
    snapshots = Dict{Int,Array{ComplexF64,3}}()
    for levelID in 2:octree.nLevels
        level = octree.levels[levelID]
        isdefined(level, :disaggG) || continue
        snapshots[levelID] = ComplexF64.(copy(level.disaggG))
    end
    return snapshots
end

function worst_level_disagg_difference(approx_snapshots, exact_snapshots)
    worst_level = 0
    worst_rel = -1.0
    worst_abs = 0.0
    worst_ref = 0.0

    for levelID in sort(collect(intersect(keys(approx_snapshots), keys(exact_snapshots))))
        approx_field = approx_snapshots[levelID]
        exact_field = exact_snapshots[levelID]
        abs_err = norm(approx_field - exact_field)
        ref_norm = norm(exact_field)
        rel_err = abs_err / (ref_norm + 1e-30)
        if rel_err > worst_rel
            worst_level = levelID
            worst_rel = rel_err
            worst_abs = abs_err
            worst_ref = ref_norm
        end
    end

    return (level = worst_level, rel = worst_rel, abs = worst_abs, ref = worst_ref)
end

function cube_disagg_difference(approx_snapshots, exact_snapshots, levelID, cubeID)
    approx_cube = @view approx_snapshots[levelID][:, :, cubeID]
    exact_cube = @view exact_snapshots[levelID][:, :, cubeID]
    abs_err = norm(approx_cube - exact_cube)
    ref_norm = norm(exact_cube)
    rel_err = abs_err / (ref_norm + 1e-30)
    return (rel = rel_err, abs = abs_err, ref = ref_norm)
end

function run_translate_downward!(octree)
    for levelID in 2:octree.nLevels
        translate!(octree.levels[levelID])
    end
    for levelID in 2:(octree.nLevels - 1)
        disaggregate_downward!(octree.levels[levelID], octree.levels[levelID + 1])
    end
end

function apply_exact_upward_chain!(octree, basis, x_phys, sorted_ids, x_range, k; quad_order = 4)
    leaf_level = octree.levels[octree.nLevels]
    leaf_exact = exact_level_agg(leaf_level, basis, x_phys, sorted_ids, x_range, k; quad_order)
    leaf_level.aggS .= leaf_exact

    exact_aggs = Dict{Int,Array{ComplexF64,3}}(octree.nLevels => ComplexF64.(copy(leaf_exact)))
    for levelID in (octree.nLevels - 1):-1:2
        parent_level = octree.levels[levelID]
        exact_parent = exact_level_agg(parent_level, basis, x_phys, sorted_ids, x_range, k; quad_order)
        parent_level.aggS = ComplexF64.(copy(exact_parent))
        exact_aggs[levelID] = ComplexF64.(copy(exact_parent))
    end

    return exact_aggs
end

function leaf_cube_receive_metrics(octree, basis, pmchw, sorted_ids, kmode, x_m, refs)
    N = num_basis(basis)
    k, eta = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)
    factor_EM = -im * k / (4π)
    factor_HM = im * k / (eta * 4π)
    ref = kmode === :k0 ? refs.k0 : refs.k1

    leaf_level = octree.levels[octree.nLevels]
    worst_total = (cube = 0, rel = -1.0, rel_e = 0.0, rel_h = 0.0, nbf = 0)

    for iCube in 1:leaf_level.nCubes
        cube = leaf_level.cubes[iCube]
        isempty(cube.bfInterval) && continue
        field = @view leaf_level.disaggG[:, :, iCube]
        r0 = cube.center
        ids = [sorted_ids[idx] for idx in cube.bfInterval]

        fast_e = zeros(ComplexF64, length(ids))
        fast_h = zeros(ComplexF64, length(ids))
        for (local_idx, bfID_sorted) in enumerate(cube.bfInterval)
            bfID = sorted_ids[bfID_sorted]
            bf = basis.functions[bfID]
            te, tm = PMCHWMLFMAOperatorModule._receive_terms(bf, basis, field, r0, k, leaf_level.poles)
            fast_e[local_idx] = tm * factor_EM
            fast_h[local_idx] = te * factor_HM
        end

        dense_e = ref.EM[ids, :] * x_m
        dense_h = ref.HM[ids, :] * x_m
        rel_e = norm(fast_e - dense_e) / (norm(dense_e) + 1e-30)
        rel_h = norm(fast_h - dense_h) / (norm(dense_h) + 1e-30)
        rel_total = norm(vcat(fast_e, fast_h) - vcat(dense_e, dense_h)) / (norm(vcat(dense_e, dense_h)) + 1e-30)

        if rel_total > worst_total.rel
            worst_total = (cube = iCube, rel = rel_total, rel_e = rel_e, rel_h = rel_h, nbf = length(ids))
        end
    end

    return worst_total
end

function receive_terms_with_rule(bf, basis, field, r0, k, poles, quad_order)
    CT = typeof(complex(one(real(k))))
    FT = real(CT)
    te = zero(CT)
    tm = zero(CT)
    JK = CT(im * k)

    poles_r̂ = [p.r̂ for p in poles.r̂sθsϕs]
    poles_θhat = [p.θhat for p in poles.r̂sθsϕs]
    poles_ϕhat = [p.ϕhat for p in poles.r̂sθsϕs]
    nPoles = length(poles_r̂)

    tri_info = get_triangles_info(basis.mesh, basis)
    gq = GaussQuadratureInfo(:Triangle, quad_order, FT)
    n_qp = length(gq.weight)

    for i_supp in 1:2
        tri_idx = bf.support[i_supp]
        tri_idx == 0 && continue

        tri = tri_info[tri_idx]
        v_all = tri.vertices
        local_edge = bf.local_edge_idx[i_supp]
        v_opp = v_all[:, local_edge]
        sign_supp = bf.signs[i_supp]

        for i_qp in 1:n_qp
            L = gq.coordinate[:, i_qp]
            r = v_all * L
            rho = r - v_opp
            r_local = r - r0
            w_f = sign_supp * bf.edge_length / 2 * gq.weight[i_qp]

            for iPole in 1:nPoles
                r̂ = poles_r̂[iPole]
                θhat = poles_θhat[iPole]
                ϕhat = poles_ϕhat[iPole]
                phase = exp(-JK * dot(r̂, r_local))

                Eθ = CT(field[iPole, 1])
                Eϕ = CT(field[iPole, 2])

                E_inc = (Eθ * θhat + Eϕ * ϕhat) * phase
                te += dot(rho, E_inc) * w_f

                rhat_cross_E = (Eθ * ϕhat - Eϕ * θhat) * phase
                tm += dot(rho, rhat_cross_E) * w_f
            end
        end
    end

    return te, tm
end

function leaf_cube_receive_with_rule(octree, basis, pmchw, sorted_ids, kmode, x_m, refs, cubeID, quad_order)
    k, eta = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)
    factor_EM = -im * k / (4π)
    factor_HM = im * k / (eta * 4π)
    ref = kmode === :k0 ? refs.k0 : refs.k1

    leaf_level = octree.levels[octree.nLevels]
    cube = leaf_level.cubes[cubeID]
    field = @view leaf_level.disaggG[:, :, cubeID]
    r0 = cube.center
    ids = [sorted_ids[idx] for idx in cube.bfInterval]

    fast_e = ComplexF64[]
    fast_h = ComplexF64[]
    for bfID_sorted in cube.bfInterval
        bfID = sorted_ids[bfID_sorted]
        bf = basis.functions[bfID]
        te, tm = receive_terms_with_rule(bf, basis, field, r0, k, leaf_level.poles, quad_order)
        push!(fast_e, tm * factor_EM)
        push!(fast_h, te * factor_HM)
    end

    dense_e = ref.EM[ids, :] * x_m
    dense_h = ref.HM[ids, :] * x_m
    rel_e = norm(fast_e - dense_e) / (norm(dense_e) + 1e-30)
    rel_h = norm(fast_h - dense_h) / (norm(dense_h) + 1e-30)
    rel_total = norm(vcat(fast_e, fast_h) - vcat(dense_e, dense_h)) / (norm(vcat(dense_e, dense_h)) + 1e-30)

    return (rel_total = rel_total, rel_e = rel_e, rel_h = rel_h, nbf = length(ids))
end

function leaf_cube_basis_detail(octree, basis, pmchw, sorted_ids, kmode, x_m, refs, cubeID)
    k, eta = kmode === :k0 ? (pmchw.k0, pmchw.eta0) : (pmchw.k1, pmchw.eta1)
    factor_EM = -im * k / (4π)
    factor_HM = im * k / (eta * 4π)
    ref = kmode === :k0 ? refs.k0 : refs.k1

    leaf_level = octree.levels[octree.nLevels]
    cube = leaf_level.cubes[cubeID]
    field = @view leaf_level.disaggG[:, :, cubeID]
    r0 = cube.center
    ids = [sorted_ids[idx] for idx in cube.bfInterval]

    rows = NamedTuple[]
    for bfID_sorted in cube.bfInterval
        bfID = sorted_ids[bfID_sorted]
        bf = basis.functions[bfID]
        te, tm = PMCHWMLFMAOperatorModule._receive_terms(bf, basis, field, r0, k, leaf_level.poles)
        fast_e = tm * factor_EM
        fast_h = te * factor_HM
        dense_e = dot(@view(ref.EM[bfID, :]), x_m)
        dense_h = dot(@view(ref.HM[bfID, :]), x_m)
        push!(rows, (
            basis = bfID,
            fast_e = fast_e,
            dense_e = dense_e,
            rel_e = abs(fast_e - dense_e) / (abs(dense_e) + 1e-30),
            fast_h = fast_h,
            dense_h = dense_h,
            rel_h = abs(fast_h - dense_h) / (abs(dense_h) + 1e-30),
        ))
    end

    return rows
end

function leaf_cube_field_difference(approx_snapshots, exact_snapshots, levelID, cubeID)
    approx_cube = @view approx_snapshots[levelID][:, :, cubeID]
    exact_cube = @view exact_snapshots[levelID][:, :, cubeID]

    rows = NamedTuple[]
    for iPole in axes(exact_cube, 1)
        approx_theta = approx_cube[iPole, 1]
        approx_phi = approx_cube[iPole, 2]
        exact_theta = exact_cube[iPole, 1]
        exact_phi = exact_cube[iPole, 2]
        abs_err = norm((approx_theta - exact_theta, approx_phi - exact_phi))
        ref_norm = norm((exact_theta, exact_phi))
        rel_err = abs_err / (ref_norm + 1e-30)
        push!(rows, (
            pole = iPole,
            abs_err = abs_err,
            rel_err = rel_err,
            approx_theta = approx_theta,
            exact_theta = exact_theta,
            approx_phi = approx_phi,
            exact_phi = exact_phi,
        ))
    end

    sort!(rows, by = row -> row.abs_err, rev = true)
    return rows
end

function analyze_case(name, budget, pmchw, basis, x_phys)
    op = PMCHWMLFMAOperator(pmchw, basis, 0.10; budget = budget)
    refs = direct_m_pass_refs(pmchw, basis, op)
    N = num_basis(basis)
    x_m = view(x_phys, (N + 1):(2N))

    println()
    @printf("=== %s ===\n", name)
    @printf("leaf_size_eff = %.5f, near_range = %d\n", op.leaf_size_eff, op.near_range)

    for kmode in (:k0, :k1)
        octree = kmode === :k0 ? op.octree0 : op.octree1
        sorted_ids = kmode === :k0 ? op.sorted_ids0 : op.sorted_ids1
        k = kmode === :k0 ? pmchw.k0 : pmchw.k1

        clear_agg!(octree)
        aggregate_leaf_pmchw!(octree, basis, x_phys, sorted_ids, (N + 1):(2N), k)

        println()
        @printf("[%s] upward exact-reintegration check (4-point source)\n", String(kmode))

        first_bad = nothing
        for levelID in (octree.nLevels - 1):-1:2
            parent_level = octree.levels[levelID]
            child_level = octree.levels[levelID + 1]
            aggregate_upward!(parent_level, child_level)
            exact_parent = exact_level_agg(parent_level, basis, x_phys, sorted_ids, (N + 1):(2N), k)
            worst = worst_cube_metrics(parent_level.aggS, exact_parent)
            @printf(
                "  level %d: worst cube = %d, rel = %.6e, abs = %.6e, ref = %.6e\n",
                levelID,
                worst.cube,
                worst.rel,
                worst.abs,
                worst.ref,
            )
            if isnothing(first_bad) && worst.rel > 1e-3
                first_bad = (level = levelID, cube = worst.cube, rel = worst.rel)
            end
        end

        if isnothing(first_bad)
            println("  first upward mismatch above 1e-3: none")
        else
            @printf(
                "  first upward mismatch above 1e-3: level %d cube %d rel %.6e\n",
                first_bad.level,
                first_bad.cube,
                first_bad.rel,
            )
        end

        run_translate_downward!(octree)
        approx_disagg_snapshots = capture_disagg_snapshots(octree)

        leaf_worst = leaf_cube_receive_metrics(octree, basis, pmchw, sorted_ids, kmode, x_m, refs)
        @printf(
            "[%s] approx leaf receive worst cube = %d, rel_total = %.6e, rel_e = %.6e, rel_h = %.6e, nbf = %d\n",
            String(kmode),
            leaf_worst.cube,
            leaf_worst.rel,
            leaf_worst.rel_e,
            leaf_worst.rel_h,
            leaf_worst.nbf,
        )

        approx_basis_rows = leaf_cube_basis_detail(octree, basis, pmchw, sorted_ids, kmode, x_m, refs, leaf_worst.cube)

        clear_agg!(octree)
        apply_exact_upward_chain!(octree, basis, x_phys, sorted_ids, (N + 1):(2N), k)
        run_translate_downward!(octree)
        exact_disagg_snapshots = capture_disagg_snapshots(octree)

        exact_upward_leaf_worst = leaf_cube_receive_metrics(octree, basis, pmchw, sorted_ids, kmode, x_m, refs)
        worst_disagg = worst_level_disagg_difference(approx_disagg_snapshots, exact_disagg_snapshots)
        cube_leaf_disagg = cube_disagg_difference(approx_disagg_snapshots, exact_disagg_snapshots, octree.nLevels, leaf_worst.cube)
        @printf(
            "[%s] exact-upward leaf receive worst cube = %d, rel_total = %.6e, rel_e = %.6e, rel_h = %.6e, nbf = %d\n",
            String(kmode),
            exact_upward_leaf_worst.cube,
            exact_upward_leaf_worst.rel,
            exact_upward_leaf_worst.rel_e,
            exact_upward_leaf_worst.rel_h,
            exact_upward_leaf_worst.nbf,
        )
        @printf(
            "[%s] approx-vs-exact-upward disaggG worst level = %d, rel = %.6e, abs = %.6e, ref = %.6e\n",
            String(kmode),
            worst_disagg.level,
            worst_disagg.rel,
            worst_disagg.abs,
            worst_disagg.ref,
        )
        @printf(
            "[%s] approx-vs-exact-upward leaf disaggG on approx worst cube %d: rel = %.6e, abs = %.6e, ref = %.6e\n",
            String(kmode),
            leaf_worst.cube,
            cube_leaf_disagg.rel,
            cube_leaf_disagg.abs,
            cube_leaf_disagg.ref,
        )

        println("  approx chain basis detail on approx worst cube:")
        for row in approx_basis_rows
            @printf(
                "    basis %d: rel_e = %.6e, rel_h = %.6e, fast_e = %.6e%+.6ei, dense_e = %.6e%+.6ei, fast_h = %.6e%+.6ei, dense_h = %.6e%+.6ei\n",
                row.basis,
                row.rel_e,
                row.rel_h,
                real(row.fast_e),
                imag(row.fast_e),
                real(row.dense_e),
                imag(row.dense_e),
                real(row.fast_h),
                imag(row.fast_h),
                real(row.dense_h),
                imag(row.dense_h),
            )
        end

        alt_rule4 = leaf_cube_receive_with_rule(octree, basis, pmchw, sorted_ids, kmode, x_m, refs, leaf_worst.cube, 4)
        alt_rule7 = leaf_cube_receive_with_rule(octree, basis, pmchw, sorted_ids, kmode, x_m, refs, leaf_worst.cube, 7)
        @printf(
            "[%s] exact-upward worst cube %d leaf receive with 4-point rule: rel_total = %.6e, rel_e = %.6e, rel_h = %.6e, nbf = %d\n",
            String(kmode),
            leaf_worst.cube,
            alt_rule4.rel_total,
            alt_rule4.rel_e,
            alt_rule4.rel_h,
            alt_rule4.nbf,
        )
        @printf(
            "[%s] exact-upward worst cube %d leaf receive with 7-point rule: rel_total = %.6e, rel_e = %.6e, rel_h = %.6e, nbf = %d\n",
            String(kmode),
            leaf_worst.cube,
            alt_rule7.rel_total,
            alt_rule7.rel_e,
            alt_rule7.rel_h,
            alt_rule7.nbf,
        )

        clear_agg!(octree)
        apply_exact_upward_chain!(octree, basis, x_phys, sorted_ids, (N + 1):(2N), k; quad_order = 3)
        run_translate_downward!(octree)
        exact3_disagg_snapshots = capture_disagg_snapshots(octree)
        exact3_leaf_worst = leaf_cube_receive_metrics(octree, basis, pmchw, sorted_ids, kmode, x_m, refs)
        exact3_cube_leaf_disagg = cube_disagg_difference(approx_disagg_snapshots, exact3_disagg_snapshots, octree.nLevels, leaf_worst.cube)
        @printf(
            "[%s] exact-upward(3pt src legacy) leaf receive worst cube = %d, rel_total = %.6e, rel_e = %.6e, rel_h = %.6e, nbf = %d\n",
            String(kmode),
            exact3_leaf_worst.cube,
            exact3_leaf_worst.rel,
            exact3_leaf_worst.rel_e,
            exact3_leaf_worst.rel_h,
            exact3_leaf_worst.nbf,
        )
        @printf(
            "[%s] approx-vs-exact-upward(3pt src legacy) leaf disaggG on approx worst cube %d: rel = %.6e, abs = %.6e, ref = %.6e\n",
            String(kmode),
            leaf_worst.cube,
            exact3_cube_leaf_disagg.rel,
            exact3_cube_leaf_disagg.abs,
            exact3_cube_leaf_disagg.ref,
        )

        field_rows = leaf_cube_field_difference(approx_disagg_snapshots, exact3_disagg_snapshots, octree.nLevels, leaf_worst.cube)
        println("  top pole differences on approx worst cube vs exact-upward(3pt src legacy):")
        for row in field_rows[1:min(5, length(field_rows))]
            @printf(
                "    pole %d: abs = %.6e, rel = %.6e, approx_theta = %.6e%+.6ei, exact_theta = %.6e%+.6ei, approx_phi = %.6e%+.6ei, exact_phi = %.6e%+.6ei\n",
                row.pole,
                row.abs_err,
                row.rel_err,
                real(row.approx_theta),
                imag(row.approx_theta),
                real(row.exact_theta),
                imag(row.exact_theta),
                real(row.approx_phi),
                imag(row.approx_phi),
                real(row.exact_phi),
                imag(row.exact_phi),
            )
        end
    end
end

pmchw, basis = make_fixture()
x_m = make_probe(num_basis(basis), :M; seed = 43)

analyze_case("default", PMCHWMLFMAErrorBudget(Float64), pmchw, basis, x_m)
analyze_case(
    "loose_near",
    PMCHWMLFMAErrorBudget(Float64; near_range_scale = 4.0, min_near_range = 4, max_near_range = 16),
    pmchw,
    basis,
    x_m,
)