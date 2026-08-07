using EMMoMSuite
using MoM_Basics
using MoM_Kernels
using LinearAlgebra

EMMoMSuite.set_frequency!(1.2e9)

mesh_file = joinpath(@__DIR__, "..", "..", "..", "deps", "fixtures", "AllinOne", "meshfiles", "plate_1dot2GHz.nas")
mesh = EMMoMSuite.read_nas_mesh(mesh_file, scale = 1.0)
swg = EMMoMSuite.SWGBasis(mesh)
perms = fill(ComplexF64(2.0 * (1 - 0.0002im)), EMMoMSuite.CoreModule.num_elements(mesh))
vefie = EMMoMSuite.VEFIE(1.2e9, perms)
ems_tetras = EMMoMSuite.get_tetrahedra_info(mesh, swg, perms)
ems_cache = EMMoMSuite.IntegralEquations.precompute_vefie_basis(vefie, ems_tetras)

MoM_Basics.setPrecision!(Float64)
MoM_Basics.SimulationParams.SHOWIMAGE = false
MoM_Kernels.inputParameters(
    frequency = 1.2e9,
    ieT = :EFIE,
    meshfilename = mesh_file,
)
MoM_Basics.updateVSBFTParams!(; sbfT = :nothing, vbfT = :SWG)

mesh_data, _ = MoM_Basics.getMeshData(mesh_file; meshUnit = :m)
_, legacy_nbf, legacy_tetras, _ = MoM_Basics.getCellsBFs(mesh_data, :SWG)
MoM_Basics.setGeosPermittivity!(legacy_tetras, ComplexF64(2.0 * (1 - 0.0002im)))

ems_source = EMMoMSuite.PlaneWave(1.2e9, 3pi / 4, pi, [-1.0, 0.0, 1.0])
legacy_source = MoM_Basics.PlaneWave(pi / 4, 0.0, 0.0, 1.0)

ems_v = EMMoMSuite.excitation_vector(vefie, ems_source, swg, perms)
legacy_v = ComplexF64.(MoM_Kernels.getExcitationVector(legacy_tetras, legacy_nbf, legacy_source))

println("ems_tetras=$(length(ems_tetras)) legacy_tetras=$(length(legacy_tetras))")
println("ems_nbf=$(EMMoMSuite.num_basis(swg)) legacy_nbf=$(legacy_nbf)")
println("center_diff_1=$(norm(ems_tetras[1].center - legacy_tetras[1].center))")
println("excitation_maxabs=$(maximum(abs.(ems_v - legacy_v))) excitation_fro=$(norm(ems_v - legacy_v))")
println("excitation_first=$(ems_v[1]) legacy_first=$(legacy_v[1])")

probe_i = zeros(ComplexF64, legacy_nbf)
probe_i[1] = 1.0 + 0.0im
theta_obs = [0.3]
phi_obs = [0.0, pi / 2]
_, ems_rcs_total, ems_rcs_db = EMMoMSuite.radarCrossSection(theta_obs, phi_obs, probe_i, swg, perms)
ems_j = EMMoMSuite.geoElectricJCal(probe_i, swg, perms)

function compute_legacy_rcs_total(coeffs, geos_info, theta_vals, phi_vals)
    j_geos = MoM_Kernels.geoElectricJCal(coeffs, geos_info, MoM_Basics.VSBFTypes.vbfType)
    rcs = zeros(Float64, length(theta_vals), length(phi_vals))
    factor = (MoM_Basics.Params.K_0 * MoM_Basics.η_0)^2 / (4 * pi)
    for (i, theta) in enumerate(theta_vals)
        theta_info = MoM_Basics.∠Info{Float64}(theta)
        for (j, phi) in enumerate(phi_vals)
            r_info = MoM_Basics.r̂θϕInfo(theta_info, MoM_Basics.∠Info{Float64}(phi))
            n_theta_phi = MoM_Kernels.raditionalIntegralNθϕCal(r_info, geos_info, j_geos)
            rcs[i, j] = factor * sum(abs2, n_theta_phi)
        end
    end
    return rcs
end

legacy_rcs_total = compute_legacy_rcs_total(probe_i, legacy_tetras, theta_obs, phi_obs)
legacy_rcs_db = 10 .* log10.(legacy_rcs_total)
legacy_j = MoM_Kernels.geoElectricJCal(probe_i, legacy_tetras, MoM_Basics.VSBFTypes.vbfType)
println("probe_current_maxabs=$(maximum(abs.(ems_j - legacy_j))) probe_current_fro=$(norm(ems_j - legacy_j))")
println("rcs_probe_total_diff=$(maximum(abs.(ems_rcs_total - legacy_rcs_total)))")
println("rcs_probe_db_0=$(ems_rcs_db[1,1]) legacy=$(legacy_rcs_db[1,1])")
println("rcs_probe_db_90=$(ems_rcs_db[1,2]) legacy=$(legacy_rcs_db[1,2])")

legacy_threshold(tet) = MoM_Basics.Params.Rsglr / sqrt(norm(tet.ε) / MoM_Basics.ε_0)

function compare_pair(label, it, js)
    legacy_t = legacy_tetras[it]
    legacy_s = legacy_tetras[js]
    ems_t = ems_tetras[it]
    ems_s = ems_tetras[js]
    cache_t = ems_cache[it]
    cache_s = ems_cache[js]

    if it == js
        legacy_block = ComplexF64.(MoM_Kernels.EFIEOnTetraSWG(legacy_t))
        ems_block = Matrix(EMMoMSuite.IntegralEquations.VEFIEModule._ordered_swg_self_kernel(vefie, ems_t))
    else
        dist = norm(legacy_t.center - legacy_s.center)
        if dist < legacy_threshold(legacy_t)
            legacy_block, _ = MoM_Kernels.EFIEOnNearTetrasSWG(legacy_t, legacy_s)
        else
            legacy_block, _ = MoM_Kernels.EFIEOnTetrasSWG(legacy_t, legacy_s)
        end
        legacy_block = ComplexF64.(legacy_block)
        ems_block = Matrix(EMMoMSuite.IntegralEquations.vefie_element_interaction_kernel(vefie, ems_t, ems_s, cache_t, cache_s))
    end

    diff = ems_block - legacy_block
    println("[$label] pair=($it,$js) maxabs=$(maximum(abs.(diff))) fro=$(norm(diff))")
    println("[$label] ems_11=$(ems_block[1,1]) legacy_11=$(legacy_block[1,1]) diff_11=$(diff[1,1])")
end

function ems_self_decomposition(tet)
    FT = Float64
    CT = ComplexF64
    k = vefie.k
    omega = 2pi * vefie.freq
    eps0 = FT(8.854187817e-12)
    div4pi = one(FT) / (4pi)
    sscg = div4pi .* EMMoMSuite.IntegralEquations.VEFIEModule.compute_SSCg(k)
    gq_tet = EMMoMSuite.Geometry.GaussQuadratureInfo(:Tetrahedron, 11, FT)
    gq_tri = EMMoMSuite.Geometry.GaussQuadratureInfo(:Triangle, 7, FT)
    rgt = tet.vertices * gq_tet.coordinate

    c1f1 = zeros(CT, 4, 4)
    cf23 = zeros(CT, 4, 4)
    f45 = zeros(CT, 4, 4)
    f6m = zeros(CT, 4, 4)

    ig_div_vs = zeros(CT, length(gq_tet.weight))
    ivecg_div_vs = Matrix{CT}(undef, 3, length(gq_tet.weight))
    for gi in eachindex(gq_tet.weight)
        rgi = @view rgt[:, gi]
        ig_v, ivecg_v = EMMoMSuite.IntegralEquations.VEFIEModule.volumeSingularityIgIvecg(rgi, tet, sscg)
        ig_div_vs[gi] = ig_v / tet.volume
        ivecg_div_vs[:, gi] .= ivecg_v / tet.volume
    end
    f3 = sum(gq_tet.weight[gi] * ig_div_vs[gi] for gi in eachindex(gq_tet.weight))

    signed_areas = tet.facesArea .* tet.bfsSign
    f4s = zeros(CT, 4)
    for ni = 1:4
        arean = signed_areas[ni]
        δκn = tet.faces[ni].δκ
        isbdn = tet.faces[ni].isbd
        if isbdn || ((δκn != 0) && (arean > 0))
            ig_s_div_s = zero(CT)
            for gi in eachindex(gq_tet.weight)
                rgi = @view rgt[:, gi]
                ig = EMMoMSuite.IntegralEquations.VEFIEModule.faceSingularityIg(
                    rgi,
                    tet.faces[ni].vertices,
                    tet.faces[ni].edgel,
                    tet.faces[ni].edgev̂,
                    tet.faces[ni].edgen̂,
                    abs(arean),
                    view(tet.facesn̂, :, ni),
                    sscg,
                )
                ig_s_div_s += gq_tet.weight[gi] * ig
            end
            f4s[ni] = ig_s_div_s / abs(arean)
        end
    end

    div_v_eps = one(FT) / (tet.volume * tet.ε * eps0)
    for ni = 1:4
        arean = signed_areas[ni]
        free_vn = @view tet.vertices[:, ni]
        δκn = tet.faces[ni].δκ
        n_global = tet.inBfsID[ni]
        for mi = 1:4
            aream = signed_areas[mi]
            free_vm = @view tet.vertices[:, mi]
            isbdm = tet.faces[mi].isbd
            aman = aream * arean
            c1 = aman * div_v_eps / (im * omega * 9)
            c3 = aman / (im * omega * eps0)

            f1 = zero(Float64)
            f2 = zero(CT)
            for gi in eachindex(gq_tet.weight)
                rgi = @view rgt[:, gi]
                ρm = rgi - free_vm
                ρmn = rgi - free_vn
                ρmρn = dot(ρm, ρmn)
                f1 += gq_tet.weight[gi] * ρmρn
                f2 += gq_tet.weight[gi] * (-dot(ρm, @view ivecg_div_vs[:, gi]) + ρmρn * ig_div_vs[gi])
            end

            c1f1[mi, ni] = c1 * f1
            cf23[mi, ni] = tet.κ * c3 * (-(k^2 / 9) * f2 + f3)
            if (δκn != 0) && (arean > 0)
                f45[mi, ni] -= δκn * c3 * f4s[ni]
            end
            if isbdm
                f45[mi, ni] -= tet.κ * c3 * f4s[mi]
            end

            if isbdm && (δκn != 0) && (arean > 0)
                f6 = if n_global == tet.inBfsID[mi]
                    EMMoMSuite.IntegralEquations.VEFIEModule._same_face_f6(tet.faces[ni], gq_tri, k, div4pi)
                else
                    rq_face_m = EMMoMSuite.IntegralEquations.VEFIEModule._triangle_points(tet.faces[mi], gq_tri)
                    rq_face_n = EMMoMSuite.IntegralEquations.VEFIEModule._triangle_points(tet.faces[ni], gq_tri)
                    acc = zero(CT)
                    for gj in eachindex(gq_tri.weight)
                        rgj = @view rq_face_m[:, gj]
                        for gi in eachindex(gq_tri.weight)
                            rgi = @view rq_face_n[:, gi]
                            acc += EMMoMSuite.IntegralEquations.fast_green_func(vefie.exp_table, norm(rgi - rgj)) * gq_tri.weight[gi] * gq_tri.weight[gj]
                        end
                    end
                    acc
                end
                f6m[mi, ni] += δκn * c3 * f6
            end
        end
    end

    return (; c1f1, cf23, f45, f6 = f6m, total = c1f1 + cf23 + f45 + f6m)
end

function legacy_self_decomposition(tet)
    FT = Float64
    CT = ComplexF64
    k = MoM_Basics.Params.K_0
    omega = 2pi * vefie.freq
    eps0 = FT(MoM_Basics.ε_0)
    tet_weights = Float64.(MoM_Basics.TetraGQInfoSglr.weight)
    tri_weights = Float64.(MoM_Basics.TriGQInfoSglr.weight)
    rgt = Float64.(MoM_Basics.getGQPTetraSglr(tet))

    c1f1 = zeros(CT, 4, 4)
    cf23 = zeros(CT, 4, 4)
    f45 = zeros(CT, 4, 4)
    f6m = zeros(CT, 4, 4)

    ig_div_vs = zeros(CT, length(tet_weights))
    ivecg_div_vs = Matrix{CT}(undef, 3, length(tet_weights))
    for gi in eachindex(tet_weights)
        rgi = @view rgt[:, gi]
        ig_v, ivecg_v = MoM_Kernels.volumeSingularityIgIvecg(rgi, tet)
        ig_div_vs[gi] = ig_v / tet.volume
        ivecg_div_vs[:, gi] .= ivecg_v / tet.volume
    end
    f3 = sum(tet_weights[gi] * ig_div_vs[gi] for gi in eachindex(tet_weights))

    f4s = zeros(CT, 4)
    for ni = 1:4
        arean = tet.facesArea[ni]
        δκn = tet.faces[ni].δκ
        isbdn = tet.faces[ni].isbd
        if isbdn || ((δκn != 0) && (arean > 0))
            ig_s_div_s = zero(CT)
            for gi in eachindex(tet_weights)
                rgi = @view rgt[:, gi]
                ig = MoM_Kernels.faceSingularityIg(rgi, tet.faces[ni], abs(arean), tet.facesn̂[:, ni])
                ig_s_div_s += tet_weights[gi] * ig
            end
            f4s[ni] = ig_s_div_s / abs(arean)
        end
    end

    div_v_eps = one(FT) / (tet.volume * tet.ε)
    for ni = 1:4
        arean = tet.facesArea[ni]
        free_vn = @view tet.vertices[:, ni]
        δκn = tet.faces[ni].δκ
        n_global = tet.inBfsID[ni]
        for mi = 1:4
            aream = tet.facesArea[mi]
            free_vm = @view tet.vertices[:, mi]
            isbdm = tet.faces[mi].isbd
            aman = aream * arean
            c1 = aman * div_v_eps / (im * omega * 9)
            c3 = aman / (im * omega * eps0 * 4pi)

            f1 = zero(Float64)
            f2 = zero(CT)
            for gi in eachindex(tet_weights)
                rgi = @view rgt[:, gi]
                ρm = rgi - free_vm
                ρmn = rgi - free_vn
                ρmρn = dot(ρm, ρmn)
                f1 += tet_weights[gi] * ρmρn
                f2 += tet_weights[gi] * (-dot(ρm, @view ivecg_div_vs[:, gi]) + ρmρn * ig_div_vs[gi])
            end

            c1f1[mi, ni] = c1 * f1
            cf23[mi, ni] = tet.κ * c3 * (-(k^2 / 9) * f2 + f3)
            if (δκn != 0) && (arean > 0)
                f45[mi, ni] -= δκn * c3 * f4s[ni]
            end
            if isbdm
                f45[mi, ni] -= tet.κ * c3 * f4s[mi]
            end

            if isbdm && (δκn != 0) && (arean > 0)
                f6 = if n_global == tet.inBfsID[mi]
                    rq_face = Float64.(MoM_Basics.getGQPTriSglr(tet.faces[ni]))
                    acc = zero(CT)
                    for gj in eachindex(tri_weights)
                        rgj = @view rq_face[:, gj]
                        for gi in eachindex(tri_weights)
                            rgi = @view rq_face[:, gi]
                            acc += MoM_Kernels.greenfunc_star(rgi, rgj) * tri_weights[gi] * tri_weights[gj]
                        end
                    end
                    acc + MoM_Kernels.singularF1(tet.faces[ni].edgel...)
                else
                    rq_face_m = Float64.(MoM_Basics.getGQPTriSglr(tet.faces[mi]))
                    rq_face_n = Float64.(MoM_Basics.getGQPTriSglr(tet.faces[ni]))
                    acc = zero(CT)
                    for gj in eachindex(tri_weights)
                        rgj = @view rq_face_m[:, gj]
                        for gi in eachindex(tri_weights)
                            rgi = @view rq_face_n[:, gi]
                            acc += MoM_Kernels.greenfunc(rgi, rgj) * tri_weights[gi] * tri_weights[gj]
                        end
                    end
                    acc
                end
                f6m[mi, ni] += δκn * c3 * f6
            end
        end
    end

    return (; c1f1, cf23, f45, f6 = f6m, total = c1f1 + cf23 + f45 + f6m)
end

function compare_self_decomposition(tet_idx)
    ems = ems_self_decomposition(ems_tetras[tet_idx])
    legacy = legacy_self_decomposition(legacy_tetras[tet_idx])
    for field in (:c1f1, :cf23, :f45, :f6, :total)
        diff = getproperty(ems, field) - getproperty(legacy, field)
        println("self_decomp[$field] tet=$tet_idx maxabs=$(maximum(abs.(diff))) fro=$(norm(diff))")
        peak = argmax(abs.(diff))
        i, j = Tuple(peak)
        println("self_decomp[$field] peak=($i,$j) ems=$(getproperty(ems, field)[i, j]) legacy=$(getproperty(legacy, field)[i, j]) diff=$(diff[i, j])")
    end
end

function compare_face_metadata(tet_idx)
    ems_tet = ems_tetras[tet_idx]
    legacy_tet = legacy_tetras[tet_idx]
    println("face_meta tet=$tet_idx")
    for face_idx = 1:4
        ems_signed_area = ems_tet.facesArea[face_idx] * ems_tet.bfsSign[face_idx]
        legacy_signed_area = legacy_tet.facesArea[face_idx]
        println("face_meta[$face_idx] area ems=$ems_signed_area legacy=$legacy_signed_area diff=$(ems_signed_area - legacy_signed_area)")
        println("face_meta[$face_idx] inBfsID ems=$(ems_tet.inBfsID[face_idx]) legacy=$(legacy_tet.inBfsID[face_idx])")
        println("face_meta[$face_idx] isbd ems=$(ems_tet.faces[face_idx].isbd) legacy=$(legacy_tet.faces[face_idx].isbd)")
        println("face_meta[$face_idx] delta_kappa ems=$(ems_tet.faces[face_idx].δκ) legacy=$(legacy_tet.faces[face_idx].δκ)")
        println("face_meta[$face_idx] vertices_diff=$(norm(ems_tet.faces[face_idx].vertices - legacy_tet.faces[face_idx].vertices))")
        println("face_meta[$face_idx] edgev_diff=$(norm(ems_tet.faces[face_idx].edgev̂ - legacy_tet.faces[face_idx].edgev̂))")
        println("face_meta[$face_idx] edgen_diff=$(norm(ems_tet.faces[face_idx].edgen̂ - legacy_tet.faces[face_idx].edgen̂))")
        println("face_meta[$face_idx] facen_diff=$(norm(ems_tet.facesn̂[:, face_idx] - legacy_tet.facesn̂[:, face_idx]))")
    end
end

function compare_hybrid_self(tet_idx)
    ems_tet = ems_tetras[tet_idx]
    legacy_block = ComplexF64.(MoM_Kernels.EFIEOnTetraSWG(legacy_tetras[tet_idx]))
    self_terms = ems_self_decomposition(ems_tet)
    mass_block = Matrix(EMMoMSuite.IntegralEquations.vefie_mass_matrix_cached(vefie, ems_tet, ems_cache[tet_idx]))
    hybrid_block = mass_block + self_terms.cf23 + self_terms.f45 + self_terms.f6
    diff = hybrid_block - legacy_block
    println("hybrid_self tet=$tet_idx maxabs=$(maximum(abs.(diff))) fro=$(norm(diff))")
    println("hybrid_self tet=$tet_idx diag ems=$(hybrid_block[4,4]) legacy=$(legacy_block[4,4]) diff=$(diff[4,4])")
end

function compare_specific_pair(it, js, label)
    legacy_far, legacy_thr, dist = legacy_use_far(legacy_tetras[it], legacy_tetras[js])
    ems_far, ems_thr, _ = ems_use_far(ems_tetras[it], ems_tetras[js])
    legacy_block = legacy_far ?
        first(MoM_Kernels.EFIEOnTetrasSWG(legacy_tetras[it], legacy_tetras[js])) :
        first(MoM_Kernels.EFIEOnNearTetrasSWG(legacy_tetras[it], legacy_tetras[js]))
    ems_block = EMMoMSuite.IntegralEquations.vefie_element_interaction_kernel(vefie, ems_tetras[it], ems_tetras[js], ems_cache[it], ems_cache[js])
    diff = Matrix(ems_block) - ComplexF64.(legacy_block)
    println("$label pair=($it,$js) dist=$dist legacy_thr=$legacy_thr ems_thr=$ems_thr legacy_far=$legacy_far ems_far=$ems_far")
    println("$label pair=($it,$js) maxabs=$(maximum(abs.(diff))) fro=$(norm(diff))")
    println("$label pair=($it,$js) diff_44=$(diff[4,4])")

    if legacy_far && ems_far
        far_alt = Matrix(
            EMMoMSuite.IntegralEquations.VEFIEModule._ordered_swg_far_kernel(
                vefie,
                ems_tetras[it],
                ems_tetras[js],
                EMMoMSuite.Geometry.GaussQuadratureInfo(:Tetrahedron, 5, Float64),
                EMMoMSuite.Geometry.GaussQuadratureInfo(:Triangle, 4, Float64),
            ),
        )
        far_alt_diff = far_alt - ComplexF64.(legacy_block)
        println("$label ordered_far_alt pair=($it,$js) maxabs=$(maximum(abs.(far_alt_diff))) fro=$(norm(far_alt_diff))")
        println("$label ordered_far_alt pair=($it,$js) diff_44=$(far_alt_diff[4,4])")
    end
end

function find_reference_pairs(ref_idx)
    near_idx = nothing
    far_idx = nothing
    ref_tet = legacy_tetras[ref_idx]
    ref_thr = legacy_threshold(ref_tet)
    for js in (ref_idx + 1):length(legacy_tetras)
        dist = norm(ref_tet.center - legacy_tetras[js].center)
        if isnothing(near_idx) && dist < ref_thr
            near_idx = js
        end
        if isnothing(far_idx) && dist > ref_thr
            far_idx = js
        end
        (!isnothing(near_idx) && !isnothing(far_idx)) && break
    end
    return ref_thr, near_idx, far_idx
end

function ems_use_far(tet_t, tet_s)
    dist = norm(tet_t.center - tet_s.center)
    threshold = EMMoMSuite.IntegralEquations.VEFIEModule._vefie_regular_threshold(vefie, tet_t)
    return dist > threshold, threshold, dist
end

function legacy_use_far(tet_t, tet_s)
    dist = norm(tet_t.center - tet_s.center)
    threshold = legacy_threshold(tet_t)
    return dist >= threshold, threshold, dist
end

function scan_reference_row(ref_idx)
    legacy_ref = legacy_tetras[ref_idx]
    ems_ref = ems_tetras[ref_idx]
    cache_ref = ems_cache[ref_idx]

    total_pairs = 0
    mismatch_pairs = 0
    max_diff = 0.0
    max_diff_pair = ref_idx
    max_mismatch_diff = 0.0
    max_mismatch_pair = ref_idx

    for js in 1:length(legacy_tetras)
        js == ref_idx && continue
        total_pairs += 1

        legacy_tet = legacy_tetras[js]
        ems_tet = ems_tetras[js]
        cache_tet = ems_cache[js]

        legacy_far, legacy_thr, dist = legacy_use_far(legacy_ref, legacy_tet)
        ems_far, ems_thr, _ = ems_use_far(ems_ref, ems_tet)

        legacy_block = if legacy_far
            first(MoM_Kernels.EFIEOnTetrasSWG(legacy_ref, legacy_tet))
        else
            first(MoM_Kernels.EFIEOnNearTetrasSWG(legacy_ref, legacy_tet))
        end
        legacy_block = ComplexF64.(legacy_block)
        ems_block = Matrix(EMMoMSuite.IntegralEquations.vefie_element_interaction_kernel(vefie, ems_ref, ems_tet, cache_ref, cache_tet))
        diff_norm = norm(ems_block - legacy_block)

        if diff_norm > max_diff
            max_diff = diff_norm
            max_diff_pair = js
        end

        if legacy_far != ems_far
            mismatch_pairs += 1
            if diff_norm > max_mismatch_diff
                max_mismatch_diff = diff_norm
                max_mismatch_pair = js
            end
        end
    end

    println("row_scan_ref=$ref_idx total_pairs=$total_pairs mismatch_pairs=$mismatch_pairs")
    println("row_scan_max_diff_pair=$max_diff_pair fro=$max_diff")
    println("row_scan_max_mismatch_pair=$max_mismatch_pair fro=$max_mismatch_diff")

    legacy_far, legacy_thr, dist = legacy_use_far(legacy_ref, legacy_tetras[max_mismatch_pair])
    ems_far, ems_thr, _ = ems_use_far(ems_ref, ems_tetras[max_mismatch_pair])
    println("row_scan_mismatch_detail pair=($ref_idx,$max_mismatch_pair) dist=$dist legacy_thr=$legacy_thr ems_thr=$ems_thr legacy_far=$legacy_far ems_far=$ems_far")
end

function assemble_ems_row(row_id)
    row = zeros(ComplexF64, legacy_nbf)
    for ti in eachindex(ems_tetras)
        tet_t = ems_tetras[ti]
        local_i = findfirst(==(row_id), tet_t.inBfsID)
        isnothing(local_i) && continue

        cache_t = ems_cache[ti]
        self_block = Matrix(EMMoMSuite.IntegralEquations.VEFIEModule._ordered_swg_self_kernel(vefie, tet_t))
        for local_j = 1:4
            col_id = tet_t.inBfsID[local_j]
            col_id == 0 && continue
            row[col_id] += self_block[local_i, local_j]
        end

        for sj in eachindex(ems_tetras)
            sj == ti && continue
            tet_s = ems_tetras[sj]
            cache_s = ems_cache[sj]
            block = Matrix(EMMoMSuite.IntegralEquations.vefie_element_interaction_kernel(vefie, tet_t, tet_s, cache_t, cache_s))
            for local_j = 1:4
                col_id = tet_s.inBfsID[local_j]
                col_id == 0 && continue
                row[col_id] += block[local_i, local_j]
            end
        end
    end
    return row
end

function assemble_legacy_row(row_id)
    row = zeros(ComplexF64, legacy_nbf)
    for ti in eachindex(legacy_tetras)
        tet_t = legacy_tetras[ti]
        local_i = findfirst(==(row_id), tet_t.inBfsID)
        isnothing(local_i) && continue

        self_block = ComplexF64.(MoM_Kernels.EFIEOnTetraSWG(tet_t))
        for local_j = 1:4
            col_id = tet_t.inBfsID[local_j]
            col_id == 0 && continue
            row[col_id] += self_block[local_i, local_j]
        end

        for sj in eachindex(legacy_tetras)
            sj == ti && continue
            tet_s = legacy_tetras[sj]
            legacy_far, _, _ = legacy_use_far(tet_t, tet_s)
            block = if legacy_far
                first(MoM_Kernels.EFIEOnTetrasSWG(tet_t, tet_s))
            else
                first(MoM_Kernels.EFIEOnNearTetrasSWG(tet_t, tet_s))
            end
            block = ComplexF64.(block)
            for local_j = 1:4
                col_id = tet_s.inBfsID[local_j]
                col_id == 0 && continue
                row[col_id] += block[local_i, local_j]
            end
        end
    end
    return row
end

function compare_global_row(row_id)
    ems_row = assemble_ems_row(row_id)
    legacy_row = assemble_legacy_row(row_id)
    diff = ems_row - legacy_row
    max_idx = argmax(abs.(diff))
    println("global_row=$row_id maxabs=$(maximum(abs.(diff))) fro=$(norm(diff))")
    println("global_row_peak_col=$max_idx ems=$(ems_row[max_idx]) legacy=$(legacy_row[max_idx]) diff=$(diff[max_idx])")
    println("global_row_diag ems=$(ems_row[row_id]) legacy=$(legacy_row[row_id]) diff=$(diff[row_id])")
    return max_idx
end

function basis_support_tetras(tetras_info, basis_id)
    support = Int[]
    for ti in eachindex(tetras_info)
        basis_id in tetras_info[ti].inBfsID && push!(support, ti)
    end
    return support
end

function trace_entry_contributions(row_id, col_id)
    row_support = basis_support_tetras(legacy_tetras, row_id)
    col_support = basis_support_tetras(legacy_tetras, col_id)
    println("trace_entry row=$row_id col=$col_id row_support=$row_support col_support=$col_support")

    total_ems = 0.0 + 0.0im
    total_legacy = 0.0 + 0.0im
    for ti in row_support
        tet_t_ems = ems_tetras[ti]
        tet_t_legacy = legacy_tetras[ti]
        local_i = findfirst(==(row_id), tet_t_legacy.inBfsID)
        cache_t = ems_cache[ti]

        for sj in col_support
            tet_s_ems = ems_tetras[sj]
            tet_s_legacy = legacy_tetras[sj]
            local_j = findfirst(==(col_id), tet_s_legacy.inBfsID)
            cache_s = ems_cache[sj]

            ems_block = if ti == sj
                Matrix(EMMoMSuite.IntegralEquations.VEFIEModule._ordered_swg_self_kernel(vefie, tet_t_ems))
            else
                Matrix(EMMoMSuite.IntegralEquations.vefie_element_interaction_kernel(vefie, tet_t_ems, tet_s_ems, cache_t, cache_s))
            end

            legacy_block = if ti == sj
                ComplexF64.(MoM_Kernels.EFIEOnTetraSWG(tet_t_legacy))
            else
                legacy_far, _, _ = legacy_use_far(tet_t_legacy, tet_s_legacy)
                block = legacy_far ? first(MoM_Kernels.EFIEOnTetrasSWG(tet_t_legacy, tet_s_legacy)) : first(MoM_Kernels.EFIEOnNearTetrasSWG(tet_t_legacy, tet_s_legacy))
                ComplexF64.(block)
            end

            ems_val = ems_block[local_i, local_j]
            legacy_val = legacy_block[local_i, local_j]
            total_ems += ems_val
            total_legacy += legacy_val
            println("trace_pair=($ti,$sj) ems=$ems_val legacy=$legacy_val diff=$(ems_val - legacy_val)")
        end
    end

    println("trace_totals ems=$total_ems legacy=$total_legacy diff=$(total_ems - total_legacy)")
end

ref_thr, near_idx, far_idx = find_reference_pairs(1)
println("legacy_threshold_ref=$ref_thr near_idx=$near_idx far_idx=$far_idx")
compare_pair("self", 1, 1)
!isnothing(near_idx) && compare_pair("near", 1, near_idx)
!isnothing(far_idx) && compare_pair("far", 1, far_idx)
compare_self_decomposition(1)
scan_reference_row(1)
peak_col = compare_global_row(1)
trace_entry_contributions(1, peak_col)

peak_support = basis_support_tetras(legacy_tetras, 1)
peak_col_support = basis_support_tetras(legacy_tetras, peak_col)
for row_tet in peak_support, col_tet in peak_col_support
    compare_specific_pair(row_tet, col_tet, "peak_pair")
end
for tet_idx in peak_support
    compare_face_metadata(tet_idx)
    compare_pair("self_peak", tet_idx, tet_idx)
    compare_self_decomposition(tet_idx)
    compare_hybrid_self(tet_idx)
end