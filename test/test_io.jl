using Test
using EMMoMSuite
using EMMoMSuite.IO
using EMMoMSuite.Geometry
using EMMoMSuite.Utilities: save_sparse_matrix, load_sparse_matrix
using HDF5
using SparseArrays
using StaticArrays

@testset "IO" begin
    # Test save_RCS_txt
    theta = [0.0, pi / 2]
    phi = [0.0, pi / 2]
    rcs_data = [1.0 2.0; 3.0 4.0]

    filename_txt = "test_rcs.txt"
    save_RCS_txt(filename_txt, theta, phi, rcs_data)

    @test isfile(filename_txt)
    rm(filename_txt)

    # Test save_results_hdf5
    filename_h5 = "test_results.h5"
    data = Dict("a" => 1, "b" => [1, 2, 3])
    save_results_hdf5(filename_h5, data)

    @test isfile(filename_h5)

    # Load back to verify
    h5open(filename_h5, "r") do file
        @test read(file["a"]) == 1
        @test read(file["b"]) == [1, 2, 3]
    end

    rm(filename_h5)

    # Test save_results_hdf5 kwargs version
    filename_h5b = "test_results_kw.h5"
    save_results_hdf5(filename_h5b; x=42, y=[7.0, 8.0])
    @test isfile(filename_h5b)
    h5open(filename_h5b, "r") do file
        @test read(file["x"]) == 42
        @test read(file["y"]) ≈ [7.0, 8.0]
    end
    rm(filename_h5b)

    # Test save_vtk (TriangleMesh)
    v1 = SVector(0.0, 0.0, 0.0)
    v2 = SVector(1.0, 0.0, 0.0)
    v3 = SVector(0.0, 1.0, 0.0)
    node = Matrix(hcat(v1, v2, v3))
    triangles = reshape([1, 2, 3], 3, 1)
    mesh = TriangleMesh(1, node, triangles)

    filename_vtk = "test_mesh" # WriteVTK adds .vtu
    outfiles = save_vtk(filename_vtk, mesh)

    @test !isempty(outfiles)
    @test isfile(outfiles[1])

    # Clean up
    for f in outfiles
        rm(f)
    end

    # Test with data (point data, same length as vertices)
    data_vec = [1.0, 2.0, 3.0]  # 3 point values
    outfiles2 = save_vtk(filename_vtk, mesh, data_vec; data_name = "PointVal")
    @test !isempty(outfiles2)
    for f in outfiles2
        rm(f)
    end

    # Test with cell data (same length as elements)
    data_cell = [99.0]  # 1 cell value
    outfiles3 = save_vtk(filename_vtk, mesh, data_cell; data_name = "CellVal")
    @test !isempty(outfiles3)
    for f in outfiles3
        rm(f)
    end

    # Test mismatched data (should warn and skip)
    data_bad = [1.0, 2.0]  # neither 3 vertices nor 1 cell
    outfiles4 = save_vtk(filename_vtk, mesh, data_bad; data_name = "BadData")
    @test !isempty(outfiles4)  # Still saves, just without data
    for f in outfiles4
        rm(f)
    end

    # Test save_vtk (TetrahedraMesh)
    nodes_tet = [
        0.0  1.0  0.0  0.0;
        0.0  0.0  1.0  0.0;
        0.0  0.0  0.0  1.0
    ]
    tetras = reshape([1, 2, 3, 4], 4, 1)
    mesh_tet = TetrahedraMesh(1, nodes_tet, tetras)

    outfiles_tet = save_vtk(filename_vtk, mesh_tet)
    @test !isempty(outfiles_tet)
    for f in outfiles_tet
        rm(f)
    end

    # TetrahedraMesh with cell data
    outfiles_tet2 = save_vtk(filename_vtk, mesh_tet, [42.0]; data_name = "TetCell")
    @test !isempty(outfiles_tet2)
    for f in outfiles_tet2
        rm(f)
    end

    # TetrahedraMesh with mismatched data (should warn)
    outfiles_tet3 = save_vtk(filename_vtk, mesh_tet, [1.0, 2.0]; data_name = "TetBad")
    @test !isempty(outfiles_tet3)
    for f in outfiles_tet3
        rm(f)
    end

    # Test HDF5Utils: save_sparse_matrix / load_sparse_matrix
    filename_sparse = "test_sparse.h5"
    # Real sparse matrix — must pass HDF5.Group, not File; create a container group
    A_real = sparse([1, 2, 3, 1], [1, 2, 3, 3], [1.0, 2.5, 3.0, 4.0], 4, 4)
    h5open(filename_sparse, "w") do file
        grp = create_group(file, "matrices")
        save_sparse_matrix(grp, "A_real", A_real)
    end
    @test isfile(filename_sparse)
    A_loaded = h5open(filename_sparse, "r") do file
        load_sparse_matrix(file["matrices"], "A_real")
    end
    @test A_loaded ≈ A_real
    @test size(A_loaded) == size(A_real)
    @test nnz(A_loaded) == nnz(A_real)
    rm(filename_sparse)

    # Complex sparse matrix
    filename_sparse2 = "test_sparse_complex.h5"
    B_complex = sparse([1, 2], [1, 2], [1.0 + 2.0im, 3.0 - 1.0im], 3, 3)
    h5open(filename_sparse2, "w") do file
        grp = create_group(file, "mats")
        save_sparse_matrix(grp, "B_complex", B_complex)
    end
    @test isfile(filename_sparse2)
    B_loaded = h5open(filename_sparse2, "r") do file
        load_sparse_matrix(file["mats"], "B_complex")
    end
    @test B_loaded ≈ B_complex
    rm(filename_sparse2)

    # load_sparse_matrix: type mismatch error (no "type" attribute)
    filename_bad = "test_sparse_bad.h5"
    h5open(filename_bad, "w") do file
        grp = create_group(file, "bad")
        g = create_group(grp, "NotAMatrix")
        g["values"] = [1, 2, 3]
        # No "type" attribute set → should throw ErrorException
    end
    @test_throws ErrorException h5open(filename_bad, "r") do file
        load_sparse_matrix(file["bad"], "NotAMatrix")
    end
    rm(filename_bad)

    # ─── Phase 17.6: VTK HexahedraMesh + complex vector + multi-field ─────────

    @testset "save_vtk HexahedraMesh" begin
        m = generate_box_volume_mesh(1.0, 1.0, 1.0, 1, 1, 1)
        @test m isa HexahedraMesh
        outfiles = save_vtk("test_hex", m)
        @test !isempty(outfiles)
        @test isfile(outfiles[1])
        for f in outfiles; rm(f); end

        # With scalar cell data
        outfiles2 = save_vtk("test_hex_cdata", m, ones(m.hexnum); data_name="vol")
        @test !isempty(outfiles2)
        for f in outfiles2; rm(f); end

        # With scalar point data
        nv = size(m.node, 2)
        outfiles3 = save_vtk("test_hex_pdata", m, ones(nv); data_name="nodal")
        @test !isempty(outfiles3)
        for f in outfiles3; rm(f); end
    end

    @testset "save_vtk complex vector field" begin
        m  = generate_sphere_mesh(1.0, 6, 12)
        nv = size(m.node, 2)

        # Complex vector at each node (point data)
        Evec = [SVector(1.0+2im, 3.0+4im, 5.0+6im) for _ in 1:nv]

        # :real_imag
        outfiles_ri = save_vtk("test_Efield_ri", m, Evec; field_name="E", save_mode=:real_imag)
        @test !isempty(outfiles_ri)
        for f in outfiles_ri; rm(f); end

        # :magnitude
        outfiles_mag = save_vtk("test_Efield_mag", m, Evec; field_name="E", save_mode=:magnitude)
        @test !isempty(outfiles_mag)
        for f in outfiles_mag; rm(f); end

        # Cell data version (one vec per triangle)
        Ecell = [SVector(1.0+0im, 0.0+0im, 0.0+0im) for _ in 1:m.trinum]
        outfiles_cell = save_vtk("test_Efield_cell", m, Ecell; field_name="Jsurf")
        @test !isempty(outfiles_cell)
        for f in outfiles_cell; rm(f); end

        # Unknown save_mode → error
        @test_throws ErrorException save_vtk("test_err", m, Evec; save_mode=:unknown_mode)
    end

    @testset "save_vtk_multi" begin
        m  = generate_box_tet_mesh(1.0, 1.0, 1.0, 2, 2, 2)
        nv = size(m.node, 2)
        nt = m.tetnum

        pd = Dict{String,Any}("temperature" => rand(nv), "pressure" => rand(nv))
        cd = Dict{String,Any}("density" => rand(nt))

        outfiles = save_vtk_multi("test_multi", m; point_data=pd, cell_data=cd)
        @test !isempty(outfiles)
        for f in outfiles; rm(f); end

        # TriangleMesh
        mt = generate_rectangle_mesh(1.0, 1.0, 3, 3)
        outfiles2 = save_vtk_multi("test_multi_tri", mt; point_data=Dict{String,Any}("x" => ones(size(mt.node,2))))
        @test !isempty(outfiles2)
        for f in outfiles2; rm(f); end

        # HexahedraMesh
        mh = generate_box_volume_mesh(1.0, 1.0, 1.0, 2, 2, 2)
        outfiles3 = save_vtk_multi("test_multi_hex", mh; cell_data=Dict{String,Any}("vol_tag"=>ones(mh.hexnum)))
        @test !isempty(outfiles3)
        for f in outfiles3; rm(f); end
    end

    @testset "save_RCS_csv" begin
        # --- 1D sweep (Nφ == 1) ---
        θs = collect(0.0:10.0:90.0) .* (π/180)        # 10 angles
        ϕs = [0.0]
        Nθ, Nϕ = length(θs), length(ϕs)
        # Synthetic RCS data: RCSθsϕs[component, iθ, iϕ]
        rcs_comp = ones(2, Nθ, Nϕ)                     # 1 m² each component
        rcs_total = rcs_comp[1,:,:] .+ rcs_comp[2,:,:]  # 2 m²
        rcs_dB    = 10 .* log10.(rcs_total)

        f1 = "test_rcs_1d.csv"
        save_RCS_csv(f1, θs, ϕs, rcs_comp, rcs_total, rcs_dB)
        @test isfile(f1)
        lines = readlines(f1)
        # header + 10 data rows
        @test length(lines) == Nθ + 1
        # parse first data row
        parts = split(lines[2], ",")
        @test length(parts) == 4   # theta_deg, rcs_theta_dBsm, rcs_phi_dBsm, rcs_total_dBsm
        @test parse(Float64, strip(parts[1])) ≈ 0.0   atol=1e-10
        @test parse(Float64, strip(parts[4])) ≈ 10*log10(2.0) atol=1e-6
        rm(f1)

        # --- 2D grid (Nφ > 1) ---
        θs2 = collect(0.0:30.0:90.0) .* (π/180)        # 4 angles
        ϕs2 = collect(0.0:90.0:360.0) .* (π/180)        # 5 angles
        Nθ2, Nϕ2 = length(θs2), length(ϕs2)
        rcs_comp2  = 2.0 .* ones(2, Nθ2, Nϕ2)
        rcs_total2 = rcs_comp2[1,:,:] .+ rcs_comp2[2,:,:]
        rcs_dB2    = 10 .* log10.(rcs_total2)

        f2 = "test_rcs_2d.csv"
        save_RCS_csv(f2, θs2, ϕs2, rcs_comp2, rcs_total2, rcs_dB2)
        @test isfile(f2)
        lines2 = readlines(f2)
        # header + Nθ2*Nϕ2 data rows
        @test length(lines2) == Nθ2 * Nϕ2 + 1
        parts2 = split(lines2[2], ",")
        @test length(parts2) == 5   # theta_deg, phi_deg, rcs_theta_dBsm, rcs_phi_dBsm, rcs_total_dBsm
        rm(f2)

        # --- linear=false (already the default, re-test explicit kwarg) ---
        f3 = "test_rcs_linear.csv"
        save_RCS_csv(f3, θs, ϕs, rcs_comp, rcs_total, rcs_dB; linear=true)
        @test isfile(f3)
        lines3 = readlines(f3)
        @test length(lines3) == Nθ + 1
        parts3 = split(lines3[2], ",")
        # linear format: theta_deg, rcs_theta_m2, rcs_phi_m2, rcs_total_m2
        @test length(parts3) == 4
        @test parse(Float64, strip(parts3[4])) ≈ 2.0   atol=1e-10
        rm(f3)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # Phase 19.2 — HexMeshIO
    # ─────────────────────────────────────────────────────────────────────────

    @testset "HexMeshIO" begin

        fixture = joinpath(@__DIR__, "fixtures", "single_hex.msh")
        @test isfile(fixture)

        # 1. read_hex_mesh returns HexahedraMesh from .msh v4 fixture
        hmesh = read_hex_mesh(fixture)
        @test hmesh isa HexahedraMesh
        @test hmesh.hexnum == 1
        @test size(hmesh.node, 2) == 8         # 8 vertices of unit cube

        # 2. Correct vertex coordinates (unit cube 0–1 in each axis)
        vmin = minimum(hmesh.node, dims=2)
        vmax = maximum(hmesh.node, dims=2)
        @test all(vmin .≈ 0.0)
        @test all(vmax .≈ 1.0)

        # 3. validate_mesh accepts the clean unit-cube mesh
        @test validate_mesh(hmesh)

        # 4. read_hex_mesh errors on a triangle mesh (not hex)
        tri_mesh_path = tempname() * ".msh"
        try
            # Write a minimal triangle .msh
            open(tri_mesh_path, "w") do io
                print(io, """\$MeshFormat
4.1 0 8
\$EndMeshFormat
\$Nodes
1 3 1 3
2 1 0 3
1
2
3
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
\$EndNodes
\$Elements
1 1 1 1
2 1 2 1
1 1 2 3
\$EndElements
""")
            end
            @test_throws ErrorException read_hex_mesh(tri_mesh_path)
        finally
            isfile(tri_mesh_path) && rm(tri_mesh_path)
        end

        # 5. unsupported format keyword throws
        @test_throws ErrorException read_hex_mesh(fixture; format=:exodus)

    end  # @testset "HexMeshIO"

end

