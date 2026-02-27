# Quick V[1] comparison
using EMSuite
using LinearAlgebra

using EMSuite.IntegralEquations.EFIEModule: EFIE
using EMSuite.IntegralEquations.EFIEModule.Impedance.Geometry: TriangleInfo, GaussQuadratureInfo, get_global_quad_points
using EMSuite.CoreModule: num_elements

const MOM_ALLINONE_DIR = joinpath(@__DIR__, "..", "..", "MoM_AllinOne")
mesh_file = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
mesh = read_nas_mesh(mesh_file, scale=1.0)
freq = 1e8
set_frequency!(freq)
basis = RWGBasis(mesh)
efie = EFIE(freq)

source = PlaneWave(freq, pi/2, pi, [0.0, 0.0, 1.0])
V_ems = excitation_vector(efie, source, basis)

bf1 = basis.functions[1]

# Compute V[1] manually using the same formula
V1_manual = 0.0im
gq3 = GaussQuadratureInfo(:Triangle, 3, Float64)
for sc in 1:2
    t_i = bf1.support[sc]
    tri_i = TriangleInfo(mesh, t_i)
    sign_val = bf1.signs[sc]
    local_edge = bf1.local_edge_idx[sc]
    area_i = tri_i.area
    
    r_quad = get_global_quad_points(tri_i, gq3)
    w_quad = gq3.weight
    
    for q in 1:length(w_quad)
        rqi = r_quad[q]
        v_free = tri_i.vertices[:, local_edge]
        rho = rqi - v_free
        f_val = sign_val * (bf1.edge_length / (2 * area_i)) * rho
        E_inc = EMSuite.CoreModule.incident_field(source, rqi)
        global V1_manual += area_i * w_quad[q] * dot(f_val, E_inc)
    end
end

println("V[1] EMSuite: $(V_ems[1])")
println("V[1] Manual:  $(V1_manual)")
println("V ratio: $(abs(V_ems[1] / V1_manual))")

# Also check V[2] and V[3]
for idx in 2:5
    bf = basis.functions[idx]
    V_manual = 0.0im
    for sc in 1:2
        t_i = bf.support[sc]
        tri_i = TriangleInfo(mesh, t_i)
        sign_val = bf.signs[sc]
        local_edge = bf.local_edge_idx[sc]
        area_i = tri_i.area
        r_quad = get_global_quad_points(tri_i, gq3)
        w_quad = gq3.weight
        for q in 1:length(w_quad)
            rqi = r_quad[q]
            v_free = tri_i.vertices[:, local_edge]
            rho = rqi - v_free
            f_val = sign_val * (bf.edge_length / (2 * area_i)) * rho
            E_inc = EMSuite.CoreModule.incident_field(source, rqi)
            V_manual += area_i * w_quad[q] * dot(f_val, E_inc)
        end
    end
    println("V[$idx] EMSuite=$(V_ems[idx]), Manual=$(V_manual), ratio=$(abs(V_ems[idx]/V_manual))")
end

println("\nSelf/Far ratio confirmed = 1.0 from previous test")
println("Near ratio = 0.0485 (EMSuite has 1.25 factor instead of 1/A_src)")
println("Near interaction fix: remove 1.25, add /A_src")
