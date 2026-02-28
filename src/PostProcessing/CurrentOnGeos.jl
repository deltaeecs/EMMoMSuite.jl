module CurrentOnGeos

using StaticArrays
using LinearAlgebra
using OffsetArrays

using ...Geometry
using ...BasisFunctions
using ...Utilities.Parameters
using ...CoreModule: num_elements, vertices, elements

export geoElectricJCal

"""
    geoElectricJCal(ICoeff, basis::RWGBasis)

Calculate weighted current on triangle patches using RWG basis functions.
J = sum(I_n * f_n)
"""
function geoElectricJCal(
    ICoeff::Vector{CT},
    basis::RWGBasis{IT,FT},
) where {IT<:Integer,FT<:Real,CT<:Complex{FT}}
    mesh = basis.mesh
    ntri = num_elements(mesh)
    Jtris = zeros(CT, 3, ntri)

    points, weights = Geometry.gaussQuadratureTri(3, FT)
    n_points = length(weights)

    verts = vertices(mesh)
    elems = elements(mesh)

    for t = 1:ntri
        Jtri = @view Jtris[:, t]

        # Get triangle vertices
        v_indices = elems[:, t]
        r1 = verts[:, v_indices[1]]
        r2 = verts[:, v_indices[2]]
        r3 = verts[:, v_indices[3]]

        # Area (needed for normalization)
        # We can compute it or use TriangleInfo if available.
        # Let's compute it.
        e1 = r2 - r1
        e3 = r1 - r3
        area = 0.5 * norm(cross(e1, -e3))

        for gi = 1:n_points
            u = points[1, gi]
            v = points[2, gi]
            w = points[3, gi]
            rgi = u * r1 + v * r2 + w * r3

            Jtritemp = zero(MVector{3,CT})

            # Loop over 3 edges of the triangle
            for k = 1:3
                bf_id = basis.basis_map[k, t]
                if bf_id == 0
                    continue
                end

                bf = basis.functions[bf_id]

                # Determine sign
                # Check if t is support[1] (+) or support[2] (-)
                if bf.support[1] == t
                    sign_val = 1.0
                    # rho = r - v_op
                    # v_op is the vertex opposite to the edge.
                    # For edge k (vertices k+1, k+2), opposite is k.
                    # Wait, local edge indexing in RWG.jl:
                    # Edge 1: v2-v3 (Opposite v1)
                    # Edge 2: v3-v1 (Opposite v2)
                    # Edge 3: v1-v2 (Opposite v3)
                    # So for edge k, opposite vertex is v_indices[k].
                    v_op = verts[:, v_indices[k]]
                    rho = rgi - v_op
                else
                    sign_val = -1.0
                    # rho = v_op - r
                    # v_op is the vertex opposite to the edge in the OTHER triangle.
                    # But f_n definition:
                    # T+: rho+ = r - v+
                    # T-: rho- = v- - r
                    # Here we are in T- (since sign is -1).
                    # So we need v- (the opposite vertex in THIS triangle).
                    v_op = verts[:, v_indices[k]]
                    rho = v_op - rgi
                end

                # f_n = l_n / (2A) * rho (with sign handled by rho definition above? No)
                # Standard def: f_n = l_n / (2A+) * rho+  (in T+)
                #               f_n = l_n / (2A-) * rho-  (in T-)
                # rho+ = r - v+
                # rho- = v- - r
                # So if we use sign_val to flip rho, we can say:
                # f_n = l_n / (2A) * (sign_val * (r - v_op)) ?
                # If T+: sign=1. rho = r - v+. Correct.
                # If T-: sign=-1. rho = -1 * (r - v-) = v- - r. Correct.

                # So: vector = rgi - v_op
                # term = I * (l / 2A) * sign * vector

                term = ICoeff[bf_id] * (bf.edge_length / (2 * area)) * sign_val * (rgi - v_op)
                Jtritemp .+= term
            end

            Jtri .+= Jtritemp * weights[gi]
        end
        # Jtri is now weighted sum.
        # But wait, raditionalIntegralNθϕCal expects Jtri to be the AVERAGE current?
        # Or the current at quadrature points?
        # raditionalIntegralNθϕCal does: JSexp .+= Jtri .* (phase * weights[gi])
        # It assumes Jtri is constant.

        # If we sum Jtritemp * weights[gi], we get integral of J over triangle (normalized by area? No).
        # int J dS = sum J(r_i) * w_i * Area.
        # Here we are summing J(r_i) * w_i.
        # So Jtri = (1/Area) * int J dS.
        # This is the average current.

        # raditionalIntegralNθϕCal multiplies by Area at the end.
        # So it computes Area * sum (J_avg * phase * w_i).
        # This is approx int J * phase dS.
        # Correct.

    end

    return Jtris
end

end
