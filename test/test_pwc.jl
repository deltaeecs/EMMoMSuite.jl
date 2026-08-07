using Test
using EMMoMSuite
using StaticArrays
using LinearAlgebra

@testset "PWC System" begin

    # ==========================================
    # Helper: Build a simple 2-tetrahedra mesh
    # ==========================================
    # Tet 1: (0,0,0)-(1,0,0)-(0,1,0)-(0,0,1)
    # Tet 2: (1,0,0)-(0,1,0)-(0,0,1)-(1,1,1)
    nodes_2tet = [
        0.0 1.0 0.0 0.0 1.0;
        0.0 0.0 1.0 0.0 1.0;
        0.0 0.0 0.0 1.0 1.0
    ]
    elements_2tet = [
        1 2;
        2 3;
        3 4;
        4 5
    ]
    tags_2tet = [1, 1]

    @testset "PWC Basis - Multi Tetrahedra" begin
        mesh = TetrahedraMesh(2, nodes_2tet, elements_2tet, tags_2tet)
        basis = PWCBasis(mesh)

        # 2 tetrahedra × 3 DOFs = 6 unknowns
        @test num_basis(basis) == 6

        # Check inBfsID mapping
        @test basis.functions[1].inBfsID == SVector(1, 2, 3)
        @test basis.functions[2].inBfsID == SVector(4, 5, 6)

        # Support mapping: global ID → tet index
        @test EMMoMSuite.CoreModule.support(basis, 1) == 1  # DOF 1 → tet 1
        @test EMMoMSuite.CoreModule.support(basis, 3) == 1  # DOF 3 → tet 1
        @test EMMoMSuite.CoreModule.support(basis, 4) == 2  # DOF 4 → tet 2
        @test EMMoMSuite.CoreModule.support(basis, 6) == 2  # DOF 6 → tet 2

        # Evaluate: unit vectors
        r = SVector(0.25, 0.25, 0.25)
        @test EMMoMSuite.CoreModule.evaluate(basis, 1, r) ≈ SVector(1.0, 0.0, 0.0)  # x̂
        @test EMMoMSuite.CoreModule.evaluate(basis, 2, r) ≈ SVector(0.0, 1.0, 0.0)  # ŷ
        @test EMMoMSuite.CoreModule.evaluate(basis, 3, r) ≈ SVector(0.0, 0.0, 1.0)  # ẑ
        @test EMMoMSuite.CoreModule.evaluate(basis, 4, r) ≈ SVector(1.0, 0.0, 0.0)  # x̂ (tet 2)
        @test EMMoMSuite.CoreModule.evaluate(basis, 5, r) ≈ SVector(0.0, 1.0, 0.0)  # ŷ (tet 2)
        @test EMMoMSuite.CoreModule.evaluate(basis, 6, r) ≈ SVector(0.0, 0.0, 1.0)  # ẑ (tet 2)

        # Volumes should be positive
        @test basis.functions[1].volume > 0
        @test basis.functions[2].volume > 0
    end

    @testset "VEFIE+PWC Assembly" begin
        mesh = TetrahedraMesh(2, nodes_2tet, elements_2tet, tags_2tet)
        basis = PWCBasis(mesh)
        nbf = num_basis(basis)  # 6

        freq = 1e9
        eps_r = 2.0 + 0.0im
        perms = fill(eps_r, num_elements(mesh))

        vefie = VEFIE(freq, perms)
        Z = assemble_impedance_matrix(vefie, basis)

        # Dimension check
        @test size(Z) == (nbf, nbf)  # (6, 6)

        # Non-zero check
        @test norm(Z) > 0

        # Diagonal blocks (self terms) should be non-zero
        # Self term for tet 1: Z[1:3, 1:3], tet 2: Z[4:6, 4:6]
        @test norm(Z[1:3, 1:3]) > 0
        @test norm(Z[4:6, 4:6]) > 0

        # Self-term diagonal should include mass matrix contribution V/(jωε)
        # Diagonal elements (xx, yy, zz of self) should have significant imaginary part
        @test abs(imag(Z[1, 1])) > 0
        @test abs(imag(Z[2, 2])) > 0

        # Off-diagonal coupling block should be non-zero (tets are touching)
        @test norm(Z[1:3, 4:6]) > 0
        @test norm(Z[4:6, 1:3]) > 0
    end

    @testset "PWC Excitation Vector" begin
        mesh = TetrahedraMesh(2, nodes_2tet, elements_2tet, tags_2tet)
        basis = PWCBasis(mesh)
        nbf = num_basis(basis)  # 6

        freq = 1e9
        source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])

        V = excitation_vector(source, basis)

        # Dimension
        @test length(V) == nbf  # 6

        # Non-zero: x-polarized plane wave should excite x-components
        # DOFs 1, 4 are x-components
        @test abs(V[1]) > 0
        @test abs(V[4]) > 0

        # y, z components should be ~0 for x-polarized wave at θ=0, ϕ=0
        @test abs(V[2]) < 1e-10 * abs(V[1])
        @test abs(V[3]) < 1e-10 * abs(V[1])
    end

    @testset "VEFIE Excitation with Operator" begin
        mesh = TetrahedraMesh(2, nodes_2tet, elements_2tet, tags_2tet)
        basis = PWCBasis(mesh)
        nbf = num_basis(basis)

        freq = 1e9
        eps_r = 2.0 + 0.0im
        perms = fill(eps_r, num_elements(mesh))
        vefie = VEFIE(freq, perms)
        source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])

        V = excitation_vector(vefie, source, basis)

        @test length(V) == nbf
        @test abs(V[1]) > 0
    end

    @testset "SCFIE+RWG+PWC Coupling" begin
        # Need mixed mesh. Try to find TriTetra.nas
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
            @warn "TriTetra.nas not found — skipping SCFIE+PWC tests"
            @test_skip true
        else
            surf_mesh, vol_mesh = read_mixed_nas_mesh(mesh_file; scale=0.001)
            basis_surf = RWGBasis(surf_mesh)
            basis_vol = PWCBasis(vol_mesh)

            n_surf = num_basis(basis_surf)
            n_vol = num_basis(basis_vol)  # 3 * num_tetrahedra
            n_total = n_surf + n_vol

            freq = 2e9
            eps_r = 2.0 * (1 - 0.001im)
            perms = fill(eps_r, num_elements(vol_mesh))

            @test n_vol == 3 * num_elements(vol_mesh)

            scfie = SCFIE(freq, perms; alpha=0.5)
            Z = assemble_impedance_matrix(scfie, basis_surf, basis_vol)

            # Full system dimensions
            @test size(Z) == (n_total, n_total)
            @test norm(Z) > 0

            # Block decomposition
            Z_SS = Z[1:n_surf, 1:n_surf]
            Z_SV = Z[1:n_surf, n_surf+1:end]
            Z_VS = Z[n_surf+1:end, 1:n_surf]
            Z_VV = Z[n_surf+1:end, n_surf+1:end]

            @test norm(Z_SS) > 0
            @test norm(Z_SV) > 0
            @test norm(Z_VS) > 0
            @test norm(Z_VV) > 0

            # Combined excitation vector
            source = PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0])
            V = excitation_vector(source, basis_surf, basis_vol)
            @test length(V) == n_total
            @test norm(V) > 0

            @testset "PWC RCS (radiation_integral_pwc)" begin
                vefie_pwc = VEFIE(freq, perms)
                Z_v = assemble_impedance_matrix(vefie_pwc, basis_vol)
                V_v = excitation_vector(vefie_pwc, PlaneWave(freq, 0.0, 0.0, [1.0, 0.0, 0.0]), basis_vol)
                I_v = Z_v \ V_v

                θs = Float64[0.0, π/2]
                ϕs = Float64[0.0]
                RCSdata, RCS_total, RCS_dB = radarCrossSection(θs, ϕs, I_v, basis_vol, perms)

                @test size(RCS_total) == (2, 1)
                @test all(isfinite, RCS_total)
                @test all(RCS_total .>= 0)
            end
        end
    end
end
