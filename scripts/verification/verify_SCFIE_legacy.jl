using EMMoMSuite
using EMMoMSuite.Geometry
using EMMoMSuite.BasisFunctions
using EMMoMSuite.IntegralEquations
using EMMoMSuite.CoreModule
using StaticArrays
using LinearAlgebra

function verify_scfie_reciprocity()
    println("==================================================")
    println("   Verification: SCFIE Reciprocity Check")
    println("==================================================")

    freq = 300e6
    c0 = 299792458.0
    omega = 2π * freq
    
    # 1. Create Surface Mesh (PEC Plate)
    # 2 Triangles
    # Nodes: (0,0,0), (1,0,0), (1,1,0), (0,1,0)
    nodes_surf = [
        0.0 1.0 1.0 0.0;
        0.0 0.0 1.0 1.0;
        0.0 0.0 0.0 0.0
    ]
    tris = [
        1 2 3; # Tri 1
        1 3 4  # Tri 2
    ]
    # Transpose to match 3xN
    tris = Matrix(tris')
    mesh_surf = TriangleMesh(2, nodes_surf, tris)
    
    # 2. Create Volume Mesh (Dielectric Tet)
    # 1 Tetrahedron
    # Nodes: (0,0,1), (1,0,1), (0,1,1), (0,0,2)
    nodes_vol = [
        0.0 1.0 0.0 0.0;
        0.0 0.0 1.0 0.0;
        1.0 1.0 1.0 2.0
    ]
    tets = [
        1 2 3 4
    ]
    tets = Matrix(tets')
    mesh_vol = TetrahedraMesh(1, nodes_vol, tets)
    
    # 3. Basis Functions
    basis_surf = RWGBasis(mesh_surf)
    basis_vol = SWGBasis(mesh_vol)
    
    n_surf = num_basis(basis_surf)
    n_vol = num_basis(basis_vol)
    println("Surface Unknowns: $n_surf")
    println("Volume Unknowns: $n_vol")
    
    # 4. Operator
    eps_r = 2.0 + 0.0im
    permittivities = [eps_r]
    scfie = SCFIE(freq, permittivities)
    
    # 5. Assembly
    println("Assembling Matrix...")
    Z = assemble_impedance_matrix(scfie, basis_surf, basis_vol)
    
    # 6. Extract Blocks
    # Z = [ Z_SS  Z_SV ]
    #     [ Z_VS  Z_VV ]
    
    Z_SV = Z[1:n_surf, n_surf+1:end]
    Z_VS = Z[n_surf+1:end, 1:n_surf]
    
    # 7. Check Reciprocity
    # Expected: Z_SV = j * omega * kappa * Z_VS^T
    # kappa = (eps - eps0) / eps = (eps_r - 1) / eps_r
    # Wait, kappa definition in VEFIE:
    # In VEFIE.jl: inv_eps = (1 - kappa) / eps0 => 1/eps = (1-kappa)/eps0 => eps0/eps = 1-kappa => kappa = 1 - eps0/eps = (eps - eps0)/eps.
    # Correct.
    
    kappa = (eps_r - 1.0) / eps_r
    factor = im * omega * kappa
    
    Z_VS_T = transpose(Z_VS)
    Z_SV_predicted = factor * Z_VS_T
    
    diff = norm(Z_SV - Z_SV_predicted)
    rel_diff = diff / norm(Z_SV)
    
    println("Z_SV norm: ", norm(Z_SV))
    println("Z_VS norm: ", norm(Z_VS))
    println("Reciprocity Factor (j*w*kappa): ", factor)
    println("Difference: ", diff)
    println("Relative Difference: ", rel_diff)
    
    if rel_diff < 1e-10
        println("SUCCESS: Reciprocity Verified!")
    else
        println("FAILURE: Reciprocity Mismatch.")
        
        # Debug first element
        println("\nDebug Element (1,1):")
        val_sv = Z_SV[1,1]
        val_vs = Z_VS[1,1] # This is Z_VS[1,1], which corresponds to Z_VS_T[1,1]
        val_pred = factor * val_vs
        println("Z_SV[1,1]: ", val_sv)
        println("Z_VS[1,1]: ", val_vs)
        println("Predicted: ", val_pred)
        println("Ratio: ", val_sv / val_vs)
    end
end

verify_scfie_reciprocity()
