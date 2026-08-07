using Test
using EMMoMSuite
using EMMoMSuite.IntegralEquations
using EMMoMSuite.PostProcessing
using EMMoMSuite.FastAlgorithms
using EMMoMSuite.BasisFunctions
using EMMoMSuite.Geometry
using LinearAlgebra

# Legacy parity tests compare against the archived MoM_Basics / MoM_Kernels
# packages. They only run when those legacy packages are available on
# LOAD_PATH (e.g. a local checkout of the legacy repos); otherwise they are
# skipped so the published package remains self-contained.
const _HAS_LEGACY_PARITY_DEPS =
    Base.find_package("MoM_Basics") !== nothing &&
    Base.find_package("MoM_Kernels") !== nothing

if _HAS_LEGACY_PARITY_DEPS
    @eval using MoM_Basics
    @eval using MoM_Kernels
else
    @info "Skipping legacy parity tests: MoM_Basics / MoM_Kernels not available on LOAD_PATH"
end

if _HAS_LEGACY_PARITY_DEPS
@testset "Legacy Parity Constants" begin
    # Constants from MoM_Basics/src/ParametersSet.jl
    # JKeta_div_16pi = 1im * k0 * eta0 / (16pi)
    # div4pi = 1 / (4pi)
    # minus_jk_div_16pi2 = -1im * k0 / (16pi^2)

    freq = 300e6
    c0 = 299792458.0
    mu0 = 4 * pi * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = 2 * pi * freq / c0
    eta = sqrt(mu0 / eps0)

    # 1) EFIE factor parity
    expected_efie_factor = im * k * eta / (16 * pi)
    efie = EFIE(freq)
    @test efie.factor ≈ expected_efie_factor atol = 1e-10

    # 2) Far-field factor parity (documented expectation)
    expected_ff_factor = -im * k * eta / (4 * pi)
    println("Expected FarField Factor: ", expected_ff_factor)

    # 3) Translation factor parity (documented expectation)
    expected_trans_factor = -im * k / (16 * pi^2)
    println("Expected Translation Factor: ", expected_trans_factor)
end

@testset "Legacy Parity SWG Tetra Face Metadata" begin
    nodes = [
        0.0 1.0 0.0 0.0 1.0
        0.0 0.0 1.0 0.0 1.0
        0.0 0.0 0.0 1.0 1.0
    ]
    elements = [
        1 2
        2 3
        3 4
        4 5
    ]
    tags = [1, 1]

    mesh = TetrahedraMesh(2, nodes, elements, tags)
    basis = SWGBasis(mesh)
    permittivities = ComplexF64[2.5 + 0im, 2.5 + 0im]
    tetras = get_tetrahedra_info(mesh, basis, permittivities)

    @test length(tetras) == 2

    internal_faces = Tuple{Int, Int}[]
    for tet_idx = 1:2
        tet = tetras[tet_idx]
        for face_idx = 1:4
            if !tet.faces[face_idx].isbd
                push!(internal_faces, (tet_idx, face_idx))
            end
        end
    end

    @test length(internal_faces) == 2

    κ = (permittivities[1] - 1) / permittivities[1]
    for tet in tetras
        for face_idx = 1:4
            signed_area = tet.facesArea[face_idx] * tet.bfsSign[face_idx]
            face = tet.faces[face_idx]
            if face.isbd
                @test face.δκ ≈ tet.bfsSign[face_idx] * κ atol = 1e-12
                @test signed_area > 0
            else
                @test face.δκ ≈ 0im atol = 1e-12
            end
        end
    end

    tet1_idx, face1_idx = internal_faces[1]
    tet2_idx, face2_idx = internal_faces[2]
    signed_area_1 = tetras[tet1_idx].facesArea[face1_idx] * tetras[tet1_idx].bfsSign[face1_idx]
    signed_area_2 = tetras[tet2_idx].facesArea[face2_idx] * tetras[tet2_idx].bfsSign[face2_idx]
    @test abs(signed_area_1) ≈ abs(signed_area_2) atol = 1e-12
    @test signed_area_1 ≈ -signed_area_2 atol = 1e-12
end

@testset "Legacy Parity VEFIE Regular Threshold" begin
    nodes = [
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
    elements = reshape([1, 2, 3, 4], 4, 1)
    tags = [1]

    mesh = TetrahedraMesh(1, nodes, elements, tags)
    basis = SWGBasis(mesh)
    permittivities = ComplexF64[2.0 * (1 - 0.0002im)]
    vefie = VEFIE(1.2e9, permittivities)
    tetras = get_tetrahedra_info(mesh, basis, permittivities)

    threshold = EMMoMSuite.IntegralEquations.VEFIEModule._vefie_regular_threshold(vefie, tetras[1])
    lambda0 = 299792458.0 / 1.2e9
    expected = 0.15 * lambda0 / sqrt(abs(permittivities[1]))

    @test threshold ≈ expected atol = 1e-12
end

@testset "Legacy Parity VEFIE Singular Kernels" begin
    EMMoMSuite.set_frequency!(1.2e9)

    mesh_file = joinpath(@__DIR__, "..", "deps", "fixtures", "AllinOne", "meshfiles", "plate_1dot2GHz.nas")
    mesh = EMMoMSuite.read_nas_mesh(mesh_file, scale = 1.0)
    basis = SWGBasis(mesh)
    permittivities = fill(ComplexF64(2.0 * (1 - 0.0002im)), EMMoMSuite.CoreModule.num_elements(mesh))
    vefie = VEFIE(1.2e9, permittivities)
    tetras = get_tetrahedra_info(mesh, basis, permittivities)
    cache = EMMoMSuite.IntegralEquations.precompute_vefie_basis(vefie, tetras)

    MoM_Basics.setPrecision!(Float64)
    MoM_Basics.SimulationParams.SHOWIMAGE = false
    MoM_Kernels.inputParameters(frequency = 1.2e9, ieT = :EFIE, meshfilename = mesh_file)
    MoM_Basics.updateVSBFTParams!(; sbfT = :nothing, vbfT = :SWG)

    mesh_data, _ = MoM_Basics.getMeshData(mesh_file; meshUnit = :m)
    _, _, legacy_tetras, _ = MoM_Basics.getCellsBFs(mesh_data, :SWG)
    MoM_Basics.setGeosPermittivity!(legacy_tetras, ComplexF64(2.0 * (1 - 0.0002im)))

    self_legacy = ComplexF64.(MoM_Kernels.EFIEOnTetraSWG(legacy_tetras[1]))
    self_ems = Matrix(EMMoMSuite.IntegralEquations.VEFIEModule._ordered_swg_self_kernel(vefie, tetras[1]))
    near_legacy = ComplexF64.(first(MoM_Kernels.EFIEOnNearTetrasSWG(legacy_tetras[1], legacy_tetras[2])))
    near_ems = Matrix(EMMoMSuite.IntegralEquations.VEFIEModule._ordered_swg_near_kernel(vefie, tetras[1], tetras[2]))
    sscg = (1 / (4 * pi)) .* EMMoMSuite.IntegralEquations.VEFIEModule.compute_SSCg(vefie.k)
    probe_point = tetras[1].vertices * EMMoMSuite.Geometry.GaussQuadratureInfo(:Tetrahedron, 11, Float64).coordinate[:, 1]
    face_legacy = MoM_Kernels.faceSingularityIg(probe_point, legacy_tetras[1].faces[1], abs(legacy_tetras[1].facesArea[1]), legacy_tetras[1].facesn̂[:, 1])
    face_ems = EMMoMSuite.IntegralEquations.VEFIEModule.faceSingularityIg(
        probe_point,
        tetras[1].faces[1].vertices,
        tetras[1].faces[1].edgel,
        tetras[1].faces[1].edgev̂,
        tetras[1].faces[1].edgen̂,
        abs(tetras[1].facesArea[1]),
        tetras[1].facesn̂[:, 1],
        sscg,
    )
    vol_legacy, ivec_legacy = MoM_Kernels.volumeSingularityIgIvecg(probe_point, legacy_tetras[1])
    vol_ems, ivec_ems = EMMoMSuite.IntegralEquations.VEFIEModule.volumeSingularityIgIvecg(probe_point, tetras[1], sscg)

    for face_idx = 1:4
        @test tetras[1].facesArea[face_idx] * tetras[1].bfsSign[face_idx] ≈ legacy_tetras[1].facesArea[face_idx] atol = 1e-12
        @test tetras[1].inBfsID[face_idx] == legacy_tetras[1].inBfsID[face_idx]
        @test tetras[1].faces[face_idx].vertices ≈ legacy_tetras[1].faces[face_idx].vertices atol = 1e-12
        @test tetras[1].faces[face_idx].edgev̂ ≈ legacy_tetras[1].faces[face_idx].edgev̂ atol = 1e-12
        @test tetras[1].faces[face_idx].edgen̂ ≈ legacy_tetras[1].faces[face_idx].edgen̂ atol = 1e-12
    end

    @test face_ems ≈ face_legacy / (4 * pi) atol = 1e-10
    @test vol_ems ≈ vol_legacy / (4 * pi) atol = 1e-10
    @test ivec_ems ≈ ivec_legacy ./ (4 * pi) atol = 1e-10
    @test norm(self_ems - self_legacy) < 1e-10
    @test norm(near_ems - near_legacy) < 1e-12
    @test norm(self_ems + Matrix(EMMoMSuite.IntegralEquations.vefie_mass_matrix_cached(vefie, tetras[1], cache[1])) - self_legacy) > norm(self_ems - self_legacy)
end
end
