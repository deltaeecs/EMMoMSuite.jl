using EMSuite
using EMSuite.Geometry
using EMSuite.BasisFunctions
using EMSuite.IntegralEquations
using EMSuite.FastAlgorithms.MLFMA
using EMSuite.FastAlgorithms.MLFMA.Aggregation
using EMSuite.FastAlgorithms.MLFMA.Translation
using EMSuite.FastAlgorithms.MLFMA.Disaggregation
using EMSuite.FastAlgorithms.MLFMA.Level
using EMSuite.FastAlgorithms.MLFMA.Octree
using LinearAlgebra
using StaticArrays

function debug_mlfma_factor()
    println("=== Debug MLFMA Factor for VEFIE ===")
    
    # 1. Setup minimal problem
    # Two tetrahedra, separated by R
    # Tet 1 at Origin.
    # Tet 2 at (0, 0, R).
    # Size L.
    
    L = 1.0 # 1 meter
    R = 10.0 # 10 meters
    
    # Vertices
    # Tet 1
    v1 = [
        0.0 1.0 0.0 0.0;
        0.0 0.0 1.0 0.0;
        0.0 0.0 0.0 1.0
    ] * L
    
    # Tet 2 (shifted by R in z)
    v2 = v1 .+ [0.0, 0.0, R]
    
    # Nodes
    nodes = hcat(v1, v2)
    # Elements
    elems = [
        [1, 2, 3, 4],
        [5, 6, 7, 8]
    ]
    
    mesh = EMSuite.Geometry.TetrahedraMesh(length(elems), nodes, hcat(elems...))
    
    # Basis
    basis = SWGBasis(mesh)
    
    # Freq
    freq = 1.5e8 # 150 MHz. lambda = 2m. 
    # R=10m = 5 lambda. Far field.
    k = 2π * freq / 299792458.0
    
    # VEFIE
    permittivities = [2.0+0im, 2.0+0im]
    vefie = VEFIE(freq, permittivities)
    
    # 2. Compute Direct Interaction Z_21 (Test 2, Source 1)
    # Tet 1 is Source (basis 1..4)
    # Tet 2 is Test (basis 5..8)
    
    # Get element infos
    tets = get_tetrahedra_info(mesh, basis, permittivities)
    tet_s = tets[1]
    tet_t = tets[2]
    
    # Compute full 4x4
    Z_ts, _ = EMSuite.IntegralEquations.VEFIEModule.vefie_element_interaction(vefie, tet_t, tet_s)
    
    println("Direct Z_ts (1,1): $(Z_ts[1,1])")
    println("Direct Norm: $(norm(Z_ts))")
    
    # 3. Compute MLFMA Interaction
    # Manual chain: Agg -> Trans -> Disagg
    
    # Setup Poles
    truncL = 5
    
    # 6-point rule
    w = 4π/6
    r̂s = [
        SVector(1.0, 0.0, 0.0), SVector(-1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0), SVector(0.0, -1.0, 0.0),
        SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -1.0)
    ]
    nPoles = 6
    
    # Compute theta hat, phi hat
    r̂sθsϕs = []
    poles_r̂ = SVector{3, Float64}[]
    poles_θhat = SVector{3, Float64}[]
    poles_ϕhat = SVector{3, Float64}[]
    Wθϕs = Float64[]
    
    for i in 1:nPoles
        r̂ = r̂s[i]
        
        # Spherical coords
        r = norm(r̂)
        θ = acos(clamp(r̂[3], -1.0, 1.0))
        ϕ = atan(r̂[2], r̂[1])
        
        # Vectors
        ct = cos(θ); st = sin(θ)
        cp = cos(ϕ); sp = sin(ϕ)
        
        # Handle gimbal lock at poles
        if abs(st) < 1e-10
            # North pole (0,0,1): theta=0. st=0.
            # Grid singularity.
            # Arbitrary check
            θhat = SVector(1.0, 0.0, 0.0)
            ϕhat = SVector(0.0, 1.0, 0.0)
        else
            θhat = SVector(ct*cp, ct*sp, -st)
            ϕhat = SVector(-sp, cp, 0.0)
        end
        
        push!(poles_r̂, r̂)
        push!(poles_θhat, θhat)
        push!(poles_ϕhat, ϕhat)
        push!(Wθϕs, w)
    end
    
    poles_struct = (r̂sθsϕs = [(r̂=poles_r̂[i], θhat=poles_θhat[i], ϕhat=poles_ϕhat[i]) for i in 1:nPoles], Wθϕs = Wθϕs)
    
    # Needs Octree info to be passed?
    # No, we pass arrays manually to internal functions.
    
    # A. Aggregation (Source Tet 1)
    # Coef = 1.0 for basis 1 (corresponding to face 1 of Tet 1)
    agg = zeros(ComplexF64, nPoles, 2, 1) # 1 cube
    
    bf_s = basis.functions[1] # Basis 1
    # Check if support is tet 1
    # bf.support is [tet_id, 0] or [0, tet_id] or [tet_id, other]
    # Here simple mesh, no shared faces?
    # Actually disjoint tets have no shared faces.
    # So each face has 1 basis function.
    # Basis 1 should be on Tet 1.
    
    # Call Agg
    center_s = [0.0, 0.0, 0.0]
    gq = EMSuite.Geometry.GaussQuadratureInfo(:Tetrahedron, 5, Float64)
    coef = 1.0 + 0im
    
    # Need access to internal function
    EMSuite.FastAlgorithms.MLFMA.Aggregation.add_radiation_pattern_swg!(agg, 1, basis, bf_s, tets, gq, k, center_s, poles_r̂, poles_θhat, poles_ϕhat, coef)
    
    radiation_pattern = agg[:, :, 1] # nPoles x 2
    
    # B. Translation
    # Center S to Center T. Vector R_vec = [0, 0, R]
    # Compute Alpha
    R_vec = [0.0, 0.0, R]
    Rab = norm(R_vec)
    x_arg = k * Rab
    
    # H2
    h2 = EMSuite.FastAlgorithms.MLFMA.Translation.spherical_h2l_array(truncL, x_arg)
    
    # Alpha
    alpha = zeros(ComplexF64, nPoles)
    const_factor = -im * k / (4 * π)
    
    for iPole in 1:nPoles
        r̂ = poles_r̂[iPole]
        cosϕ = clamp(dot(r̂, R_vec/Rab), -1.0, 1.0)
        Pl = EMSuite.FastAlgorithms.MLFMA.Translation.collectPl(truncL, cosϕ)
        
        val = 0.0 + 0im
        j_term = im
        for l in 0:truncL
             j_term *= -im
             val += j_term * (2*l+1) * h2[l+1] * Pl[l+1]
        end
        alpha[iPole] = val * const_factor * poles_struct.Wθϕs[iPole]
    end
    
    # Received Field Spectrum
    # incoming = alpha .* radiation_pattern
    incoming = zeros(ComplexF64, nPoles, 2)
    for pol in 1:2
        incoming[:, pol] = alpha .* radiation_pattern[:, pol]
    end
    
    # C. Disaggregation (Test Tet 2)
    bf_t = basis.functions[5]
    center_t = [0.0, 0.0, R]
    Z_mlfma = zeros(ComplexF64, length(basis.functions))
    EMSuite.FastAlgorithms.MLFMA.Disaggregation.add_received_field_swg!(Z_mlfma, 5, basis, bf_t, tets, gq, vefie, k, vefie.eta, center_t, poles_r̂, poles_θhat, poles_ϕhat, incoming)
    
    val_mlfma = Z_mlfma[5]
    println("MLFMA Z_ts (1,1): $val_mlfma")
    
    # 4. Compare
    ratio = abs(val_mlfma) / abs(Z_ts[1,1])
    println("Ratio: $ratio")
    
end

debug_mlfma_factor()
