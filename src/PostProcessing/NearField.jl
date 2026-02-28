module NearField

using StaticArrays
using LinearAlgebra
using ...CoreModule
using ...Geometry
using ...BasisFunctions
using ...IntegralEquations.Kernels
using ...Utilities.Parameters

export calculate_near_field

"""
    calculate_near_field(points, basis, I_coeffs)

Calculate the Electric Field at observation points.
E(r) = -j*omega*mu * A(r) - grad(Phi(r))
     = -j*omega*mu * sum(I_n * int(f_n * G)) + 1/(j*omega*eps) * sum(I_n * int(div(f_n) * grad(G)))

Arguments:
- `points`: Vector of observation points (SVector{3, FT}).
- `basis`: RWG basis.
- `I_coeffs`: Current coefficients.
"""
function calculate_near_field(
    points::Vector{SVector{3,FT}},
    basis::RWGBasis{IT,FT},
    I_coeffs::Vector{Complex{FT}},
) where {IT,FT}
    num_points = length(points)
    E_field = zeros(SVector{3,Complex{FT}}, num_points)

    omega = get_omega()
    c0 = 299792458.0
    mu0 = 4π * 1e-7
    eps0 = 1.0 / (c0^2 * mu0)
    k = omega / c0

    # Precompute constants
    const_A = -im * omega * mu0
    const_Phi = 1.0 / (im * omega * eps0)

    # Quadrature for triangles
    gq = GaussQuadratureInfo(:Triangle, 3, FT)

    mesh = basis.mesh

    # Loop over basis functions
    for n = 1:num_basis(basis)
        In = I_coeffs[n]
        if abs(In) < 1e-12
            continue
        end

        bf = basis.functions[n]

        # Loop over support triangles
        for k_supp = 1:2
            tri_idx = bf.support[k_supp]
            if tri_idx == 0
                continue
            end

            sign = bf.signs[k_supp]

            # Get triangle info
            tri = get_triangle_info(mesh, basis, tri_idx)

            # Divergence of f_n is constant on triangle
            # div(f_n) = +/- l_n / A_n
            div_f = sign * bf.edge_length / tri.area

            # Quadrature points on triangle
            q_points = get_global_quad_points(tri, gq)
            q_weights = gq.weight .* tri.area

            # Contribution to A and Phi
            for (qp, w) in zip(q_points, q_weights)
                # Value of f_n at qp
                # f_n = sign * l_n / (2*A) * rho
                # rho = r - v_op
                # We need v_op (opposite vertex)
                # In get_triangle_info, we don't explicitly identify v_op for this basis
                # But we can deduce it.
                # bf.local_edge_idx[k] is the local index of the edge.
                # The opposite vertex is the one with that local index.
                local_op_idx = bf.local_edge_idx[k_supp]
                v_op = tri.vertices[:, local_op_idx]

                rho = qp - v_op
                f_val = (sign * bf.edge_length / (2 * tri.area)) * rho

                for i = 1:num_points
                    obs = points[i]
                    G = green_function_free_space(obs, qp, k)
                    grad_G = grad_green_function_free_space(obs, qp, k)

                    # E = -j*w*mu*A - grad(Phi)
                    # Phi = (1/eps) * int(rho * G)
                    # rho = (-1/j*w) * div(J)
                    # Phi = (-1/j*w*eps) * int(div(J) * G)
                    # grad(Phi) = (-1/j*w*eps) * int(div(J) * grad(G))
                    # -grad(Phi) = (1/j*w*eps) * int(div(J) * grad(G))
                    # const_Phi = 1/(j*w*eps)
                    # So -grad(Phi) = const_Phi * div_f * grad_G

                    term_A = const_A * f_val * G
                    term_Phi = const_Phi * div_f * grad_G

                    E_field[i] += In * (term_A + term_Phi) * w
                end
            end
        end
    end

    return E_field
end

end
