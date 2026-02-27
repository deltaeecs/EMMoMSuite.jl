using Test
using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.IntegralEquations.SCFIEModule: scfie_coupling_interaction
using EMSuite.FastAlgorithms.MLFMA
using LinearAlgebra
using SparseArrays

@testset "SCFIE System" begin
    # Locate mesh file
    mesh_file = ""
    for candidate in [
        joinpath(@__DIR__, "..", "..", "MoM_Basics", "meshfiles", "TriTetra.nas"),
        joinpath(@__DIR__, "..", "..", "MoM_Kernels", "meshfiles", "TriTetra.nas"),
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
            @test rel_err < 1e-12  # Machine precision for pure near-field
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
    end
end
