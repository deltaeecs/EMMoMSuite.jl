using Test
using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.IntegralEquations.SCFIEModule: scfie_coupling_interaction, scfie_sv_only_interaction
using EMMoMSuite.FastAlgorithms.MLFMA
using LinearAlgebra
using SparseArrays

@testset "SCFIE System" begin
    # Locate mesh file
    mesh_file = ""
    for candidate in [
        joinpath(@__DIR__, "..", "deps", "fixtures", "Basics", "meshfiles", "TriTetra.nas"),
        joinpath(@__DIR__, "..", "deps", "fixtures", "Kernels", "meshfiles", "TriTetra.nas"),
    ]
        if isfile(candidate)
            mesh_file = candidate
            break
        end
    end
    
    if isempty(mesh_file)
        @warn "TriTetra.nas not found — skipping SCFIE tests"
        @test_skip true
    else
        surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001)
        basis_surf = RWGBasis(surf_mesh)
        basis_vol = SWGBasis(vol_mesh)
        n_surf = num_basis(basis_surf)
        n_vol = num_basis(basis_vol)
        n_total = n_surf + n_vol
        
        freq = 2e9
        eps_r = 2.0 * (1 - 0.001im)
        perms = fill(eps_r, num_elements(vol_mesh))
        
        @testset "Coupling reciprocity" begin
            # Verify Z_SV/κ ≈ Z_VS^T (reciprocity of L-operator)
            scfie = SCFIE(freq, perms; alpha=0.5)
            tris = get_triangles_info(surf_mesh, basis_surf)
            tetras = get_tetrahedra_info(vol_mesh, basis_vol, perms)
            
            tri = tris[1]
            tet = tetras[1]
            Z_sv, Z_vs = scfie_coupling_interaction(scfie, tri, tet)
            
            κ = tet.κ
            Z_sv_norm = Z_sv ./ κ
            recip_err = norm(Z_sv_norm - transpose(Z_vs)) / norm(Z_sv_norm)
            @test recip_err < 1e-14
        end
        
        @testset "Z_VS kappa_inv derivation" begin
            # I-7 resolution: scfie_sv_only_interaction returns Z_sv with κ embedded
            # (c1_sv = im*ω*μ₀*κ). To obtain Z_VS, divide by κ:
            #   Z_vs[n,m] = Z_sv[m,n] / κ   (reciprocity identity)
            # This MUST equal the Z_vs block from scfie_coupling_interaction.
            scfie = SCFIE(freq, perms; alpha=0.5)
            tris   = get_triangles_info(surf_mesh, basis_surf)
            tetras = get_tetrahedra_info(vol_mesh, basis_vol, perms)

            tri = tris[1]
            tet = tetras[1]

            # Full interaction (reference)
            Z_sv_full, Z_vs_full = scfie_coupling_interaction(scfie, tri, tet)

            # Optimized path: only Z_sv, then derive Z_vs via kappa_inv
            Z_sv_opt = scfie_sv_only_interaction(scfie, tri, tet)
            κ      = tet.κ
            κ_inv  = iszero(κ) ? zero(eltype(Z_sv_opt)) : one(eltype(Z_sv_opt)) / κ
            # Z_vs[n,m] = Z_sv[m,n] / κ  →  Z_vs_derived[n,m] = Z_sv_opt[m,n] * κ_inv
            Z_vs_derived = transpose(Z_sv_opt) .* κ_inv  # 4×3

            # Z_sv must match exactly
            @test norm(Z_sv_opt - Z_sv_full) / norm(Z_sv_full) < 1e-14

            # Z_vs derived must match full Z_vs
            @test norm(Z_vs_derived - Z_vs_full) / norm(Z_vs_full) < 1e-14
        end

        @testset "SCFIE Direct assembly" begin
            scfie = SCFIE(freq, perms; alpha=0.5)
            Z = assemble_impedance_matrix(scfie, basis_surf, basis_vol)
            
            @test size(Z) == (n_total, n_total)
            @test norm(Z) > 0
            
            # Block structure
            Z_SS = Z[1:n_surf, 1:n_surf]
            Z_SV = Z[1:n_surf, n_surf+1:end]
            Z_VS = Z[n_surf+1:end, 1:n_surf]
            Z_VV = Z[n_surf+1:end, n_surf+1:end]
            
            @test norm(Z_SS) > 0
            @test norm(Z_SV) > 0
            @test norm(Z_VS) > 0
            @test norm(Z_VV) > 0
        end
        
        @testset "SCFIE MLFMA near-field" begin
            scfie = SCFIE(freq, perms; alpha=0.5)
            Z_direct = assemble_impedance_matrix(scfie, basis_surf, basis_vol)
            
            # At 2GHz, 0.1m mesh is < 1λ → 100% near-field
            lambda = 299792458.0 / freq
            leaf_size = 0.25 * lambda
            mlfma_op = MLFMAOperator(scfie, [basis_surf, basis_vol], leaf_size)
            
            x = randn(ComplexF64, n_total)
            y_direct = Z_direct * x
            y_mlfma = mlfma_op * x
            
            rel_err = norm(y_direct - y_mlfma) / norm(y_direct)
            # MLFMAOperator 近场路径经 FastExp 表格插值，精度受插值截断限制（~1e-6），
            # 因此近场门限取 1e-4（4 位有效数字），而非机器精度。
            @test rel_err < 1e-4
        end
        
        @testset "SCFIE MLFMA far-field" begin
            # 4GHz → λ=0.075m, mesh=0.1m≈1.33λ → some far-field interactions
            freq2 = 4e9
            perms2 = fill(eps_r, num_elements(vol_mesh))
            scfie2 = SCFIE(freq2, perms2; alpha=0.5)
            
            Z_direct2 = assemble_impedance_matrix(scfie2, basis_surf, basis_vol)
            
            lambda2 = 299792458.0 / freq2
            leaf_size2 = 0.25 * lambda2
            mlfma_op2 = MLFMAOperator(scfie2, [basis_surf, basis_vol], leaf_size2)
            
            x2 = randn(ComplexF64, n_total)
            y_d2 = Z_direct2 * x2
            y_m2 = mlfma_op2 * x2
            
            rel_err2 = norm(y_d2 - y_m2) / norm(y_d2)
            @test rel_err2 < 0.02  # 2% tolerance for far-field
        end

        @testset "SWG RCS (radiation_integral_swg)" begin
            scfie = SCFIE(freq, perms; alpha=0.5)
            Z = assemble_impedance_matrix(scfie, basis_surf, basis_vol)
            source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
            V = excitation_vector(source, basis_surf, basis_vol)
            I_coeffs = Z \ V
            I_vol = I_coeffs[n_surf+1:end]

            θs = Float64[0.0, π/2]
            ϕs = Float64[0.0]
            RCSdata, RCS_total, RCS_dB = radarCrossSection(θs, ϕs, I_vol, basis_vol, perms)

            @test size(RCS_total) == (2, 1)
            @test all(isfinite, RCS_total)
            @test all(RCS_total .>= 0)
            @test size(RCS_dB) == (2, 1)
        end
    end
end
