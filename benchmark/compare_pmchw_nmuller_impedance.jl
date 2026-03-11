using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))

using EMSuite
using LinearAlgebra
using Printf

function feed_current_from_half(I_2N, basis, feed_edges, offset)
    N = num_basis(basis)
    current = zero(ComplexF64)
    for idx in feed_edges
        1 <= idx <= N || continue
        current += I_2N[offset + idx] * basis.functions[idx].edge_length
    end
    return current
end

function safe_impedance(current)
    return iszero(current) ? ComplexF64(Inf, Inf) : (1.0 + 0im) / current
end

function unit_edge_rhs(basis, feed_edges; offset = 0)
    N = num_basis(basis)
    V = zeros(ComplexF64, 2N)
    for idx in feed_edges
        1 <= idx <= N || continue
        V[offset + idx] = basis.functions[idx].edge_length
    end
    return V
end

function evaluate_nmuller_case(label, Z_nmuller, V_case, basis, feed_edges)
    I_case = Z_nmuller \ V_case
    front_current = feed_current_from_half(I_case, basis, feed_edges, 0)
    back_current = feed_current_from_half(I_case, basis, feed_edges, num_basis(basis))
    return (
        label = label,
        front_current = front_current,
        back_current = back_current,
        Zin_front = safe_impedance(front_current),
        Zin_back = safe_impedance(back_current),
    )
end

function print_case(result)
    @printf("  case = %s\n", result.label)
    @printf("    I_front = %+.6e + j(%+.6e)\n", real(result.front_current), imag(result.front_current))
    @printf("    I_back  = %+.6e + j(%+.6e)\n", real(result.back_current), imag(result.back_current))
    @printf("    Zin(front) = %+.6e + j(%+.6e)\n", real(result.Zin_front), imag(result.Zin_front))
    @printf("    Zin(back)  = %+.6e + j(%+.6e)\n", real(result.Zin_back), imag(result.Zin_back))
end

function parse_cli_config(args)
    preset = isempty(args) ? "small" : lowercase(args[1])
    if preset == "small"
        return (
            label = "small",
            freq = 300e6,
            eps_r = 4.0,
            mu_r = 1.0,
            radius = 0.1,
            n_theta = 4,
            n_phi = 6,
            feed_edges = [1],
        )
    elseif preset == "medium"
        return (
            label = "medium",
            freq = 300e6,
            eps_r = 4.0,
            mu_r = 1.0,
            radius = 0.1,
            n_theta = 6,
            n_phi = 10,
            feed_edges = [1],
        )
    end

    error("unsupported preset: $preset (expected: small or medium)")
end

function run_comparison(; label, freq, eps_r, mu_r, radius, n_theta, n_phi, feed_edges)
    mesh = generate_sphere_mesh(radius, n_theta, n_phi)
    basis = RWGBasis(mesh)
    pmchw = PMCHW(freq, eps_r, mu_r)
    nmuller = NMuller(freq, eps_r, mu_r)
    feed = DeltaGapSource(freq, feed_edges, 1.0 + 0im)

    Z_pmchw = assemble_impedance_matrix(pmchw, basis)
    Z_nmuller = assemble_impedance_matrix(nmuller, basis)
    V_pmchw = excitation_vector(pmchw, feed, basis)
    V_nmuller = excitation_vector(nmuller, feed, basis)

    I_pmchw = Z_pmchw \ V_pmchw
    I_nmuller = Z_nmuller \ V_nmuller
    Zin_pmchw = input_impedance(pmchw, feed, I_pmchw, basis)
    Zin_nmuller = input_impedance(nmuller, feed, I_nmuller, basis)

    feed_current_front = feed_current_from_half(I_nmuller, basis, feed_edges, 0)
    feed_current_back = feed_current_from_half(I_nmuller, basis, feed_edges, num_basis(basis))
    Zin_nmuller_front = safe_impedance(feed_current_front)
    Zin_nmuller_back = safe_impedance(feed_current_back)

    diag_results = [
        evaluate_nmuller_case("implemented heuristic", Z_nmuller, V_nmuller, basis, feed_edges),
        evaluate_nmuller_case("unit front-half only", Z_nmuller, unit_edge_rhs(basis, feed_edges; offset = 0), basis, feed_edges),
        evaluate_nmuller_case("unit back-half only", Z_nmuller, unit_edge_rhs(basis, feed_edges; offset = num_basis(basis)), basis, feed_edges),
    ]

    rel_gap = abs(Zin_nmuller - Zin_pmchw) / (abs(Zin_pmchw) + 1e-30)
    mag_ratio = abs(Zin_nmuller) / (abs(Zin_pmchw) + 1e-30)

    println("PMCHW vs N-Muller impedance comparison")
    @printf("  preset = %s\n", label)
    @printf("  N = %d (2N = %d)\n", num_basis(basis), 2 * num_basis(basis))
    @printf("  freq = %.3f MHz, eps_r = %.3f, radius = %.3f m\n", freq / 1e6, real(eps_r), radius)
    @printf("  Zin(PMCHW)   = %+.6e + j(%+.6e)\n", real(Zin_pmchw), imag(Zin_pmchw))
    @printf("  Zin(NMuller) = %+.6e + j(%+.6e)\n", real(Zin_nmuller), imag(Zin_nmuller))
    @printf("  Zin(NMuller, front current) = %+.6e + j(%+.6e)\n", real(Zin_nmuller_front), imag(Zin_nmuller_front))
    @printf("  Zin(NMuller, back current)  = %+.6e + j(%+.6e)\n", real(Zin_nmuller_back), imag(Zin_nmuller_back))
    @printf("  rel gap      = %.6e\n", rel_gap)
    @printf("  |Z| ratio     = %.6e\n", mag_ratio)
    println("  diagnostic cases:")
    for result in diag_results
        print_case(result)
    end
end

config = parse_cli_config(ARGS)
run_comparison(; config...)