using EMSuite
using EMSuite.FastAlgorithms.MLFMA.Interpolation
using EMSuite.FastAlgorithms.MLFMA.OctreeBuilder

function check_weights()
    L = 10
    # We need to call levelIntegralInfoCal or similar to get poles
    # But levelIntegralInfoCal takes cubeEdgel.
    # Let's just use the internal functions if possible, or call levelIntegralInfoCal with dummy.
    
    # Reverse engineer L
    # L = 2.72 + ...
    # Let's just pick a cubeEdgel that gives L=10.
    # Or just call octreeXWNCal directly.
    
    # From Interpolation.jl:
    # levelIntegralInfoCal calls:
    # Xs, Ws = octreeXWNCal(0.0, pi, L, :glq) (Theta)
    # Phis, Wphis = octreeXWNCal(0.0, 2pi, 2*L, :uni) (Phi)
    
    Xs, Ws = octreeXWNCal(0.0, pi, L, :glq)
    Phis, Wphis = octreeXWNCal(0.0, 2*pi, 2*L, :uni)
    
    total_weight = 0.0
    for i in 1:length(Ws)
        theta = Xs[i]
        w_theta = Ws[i]
        # Jacobian sin(theta) is included in Ws?
        # No, octreeXWNCal returns weights for the interval.
        # For GLQ on [0, pi], the weight is for integral f(x) dx.
        # Sphere integral is int f(theta, phi) sin(theta) dtheta dphi.
        # So we need to multiply by sin(theta).
        
        for j in 1:length(Wphis)
            w_phi = Wphis[j]
            total_weight += w_theta * w_phi * sin(theta)
        end
    end
    
    println("Total Weight: $total_weight")
    println("Expected: $(4*pi)")
end

check_weights()
