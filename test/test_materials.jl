using Test
using EMMoMSuite
using EMMoMSuite.CoreModule.Materials
using EMMoMSuite.CoreModule.Constants

@testset "Materials" begin
    @testset "PEC" begin
        m = PEC()
        @test permittivity(m) == Inf
        @test permeability(m) == Constants.mu0
    end

    @testset "Dielectric" begin
        # Vacuum
        d0 = Dielectric(1.0, 1.0)
        @test permittivity(d0) ≈ Constants.eps0
        @test permeability(d0) ≈ Constants.mu0
        @test impedance(d0) ≈ Constants.eta0

        # Custom material
        er = 4.0
        ur = 2.0
        d = Dielectric(er, ur)
        @test permittivity(d) ≈ er * Constants.eps0
        @test permeability(d) ≈ ur * Constants.mu0
        @test impedance(d) ≈ sqrt(ur/er) * Constants.eta0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 19.3 — MaterialLibrary
# ─────────────────────────────────────────────────────────────────────────────

@testset "MaterialLibrary" begin

    # 1. Isotropic construction (real and complex εr)
    iso = Isotropic(4.4 * (1 - 0.02im), 1.0)
    @test real(iso.εr) ≈ 4.4
    @test imag(iso.εr) ≈ -0.088
    @test iso.μr ≈ complex(1.0)

    # 2. Anisotropic with identity tensors
    using LinearAlgebra
    using StaticArrays
    eps_id = Matrix{Float64}(I, 3, 3)
    mu_id  = Matrix{Float64}(I, 3, 3)
    aniso  = Anisotropic(eps_id, mu_id)
    @test aniso.ε_tensor ≈ SMatrix{3,3,ComplexF64}(1.0I)
    @test aniso.μ_tensor ≈ SMatrix{3,3,ComplexF64}(1.0I)

    # 3. MaterialEntry construction
    entry = MaterialEntry("test_mat", Isotropic(2.1), (1e6, 20e9), v"1.0.0", "model")
    @test entry.name            == "test_mat"
    @test entry.version         == v"1.0.0"
    @test entry.frequency_range == (1e6, 20e9)

    # 4. MaterialLibrary add / get (in-range)
    lib = MaterialLibrary()
    @test isempty(lib.entries)
    add_material!(lib, entry)
    @test haskey(lib.entries, "test_mat")
    m = get_material(lib, "test_mat", 5e9)   # 5 GHz — in range
    @test m isa Isotropic
    @test real(m.εr) ≈ 2.1

    # 5. get_material: out-of-range frequency throws DomainError
    @test_throws DomainError get_material(lib, "test_mat", 200e9)   # 200 GHz > 20 GHz
    @test_throws DomainError get_material(lib, "test_mat", 1e3)     # 1 kHz < 1 MHz

    # 6. get_material: unknown material throws KeyError
    @test_throws KeyError get_material(lib, "no_such_material", 1e9)

    # 7. Built-in library: 5 standard materials present
    blib = load_builtin_library()
    @test haskey(blib.entries, "pec")
    @test haskey(blib.entries, "vacuum")
    @test haskey(blib.entries, "fr4_standard")
    @test haskey(blib.entries, "rogers_ro4003c")
    @test haskey(blib.entries, "silicon")

    pec_m = get_material(blib, "pec", 10e9)
    @test pec_m isa Isotropic
    @test isinf(real(pec_m.εr))

    vac_m = get_material(blib, "vacuum", 1e9)
    @test real(vac_m.εr) ≈ 1.0
    @test real(vac_m.μr) ≈ 1.0

    fr4_m = get_material(blib, "fr4_standard", 5e9)
    @test real(fr4_m.εr) ≈ 4.4 rtol=0.02
    @test imag(fr4_m.εr) < 0.0   # lossy: negative imag (e^{-jωt} CEM convention)

    # 8. JLD2 round-trip: save → load → same data
    tmp = tempname() * ".jld2"
    try
        save_library(blib, tmp)
        @test isfile(tmp)
        blib2 = load_library(tmp)
        @test length(blib2.entries) == length(blib.entries)
        @test haskey(blib2.entries, "fr4_standard")
        fr4_rt = get_material(blib2, "fr4_standard", 5e9)
        @test fr4_rt.εr ≈ fr4_m.εr
    finally
        isfile(tmp) && rm(tmp)
    end

end  # @testset "MaterialLibrary"

# ─────────────────────────────────────────────────────────────────────────────
# Phase 19.4 — DispersionModels
# ─────────────────────────────────────────────────────────────────────────────

@testset "DispersionModels" begin

    # 1. DebyeModel — water at 2.45 GHz  (one-pole, 25°C)
    #    Params: ε_s = 78.38, ε_inf = 4.9, τ = 8.1 ps
    #    Expected: εr ≈ 77 - j9  (validated against published Debye fits)
    water = DebyeModel(4.9, [DebyePole(78.38 - 4.9, 8.1e-12)])
    ω_2p45 = 2π * 2.45e9
    ε_w = eval_permittivity(water, ω_2p45)
    @test real(ε_w) ≈ 77.0 atol=3.0
    @test imag(ε_w) ≈ -9.0 atol=3.0   # negative: loss under e^{-jωt}

    # 2. DrudeModel — Au at 1 THz: large negative real εr
    #    Params (Johnson & Christy): ωp=1.37e16, γ=4.07e13, ε_inf=9.5
    gold = DrudeModel(9.5, 1.37e16, 4.07e13)
    ω_1thz = 2π * 1e12
    ε_au = eval_permittivity(gold, ω_1thz)
    @test real(ε_au) < -1e4          # large negative real (metallic)
    @test abs(imag(ε_au)) > 1e3      # significant imaginary (lossy)

    # 3. LorentzModel — construction + scalar eval returns ComplexF64
    lor = LorentzModel(1.0, [LorentzPole(2.0, 2π * 5e9, 2π * 0.1e9)])
    ε_lor = eval_permittivity(lor, 2π * 4.9e9)   # near resonance
    @test ε_lor isa ComplexF64
    @test real(ε_lor) > 1.0   # real part enhanced near resonance

    # 4. Vectorised eval matches scalar loop (correctness)
    ωs          = collect(range(2π * 1e9, 2π * 10e9, length=50))
    scalar_vals = [eval_permittivity(water, ω) for ω in ωs]
    vec_vals    = eval_permittivity(water, ωs)
    @test vec_vals ≈ scalar_vals

    # 5. Multi-pole DebyeModel: two identical poles double contribution above ε_inf
    water2 = DebyeModel(water.ε_inf, [water.poles[1], water.poles[1]])
    ε2 = eval_permittivity(water2, ω_2p45)
    ε1 = eval_permittivity(water,  ω_2p45)
    @test real(ε2) - water2.ε_inf ≈ 2 * (real(ε1) - water.ε_inf) rtol=1e-10

end  # @testset "DispersionModels"

