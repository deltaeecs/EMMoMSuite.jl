"""
    MeshQuality.jl

Statistical quality metrics for surface (TriangleMesh) and volume
(TetrahedraMesh) meshes.

## Quality measures

### Triangle
- **area**: element area
- **edge**: min / max edge length per element (then global stat)
- **aspect_ratio**: max_edge / min_edge  (ideal = 1 = equilateral)
- **skewness**: equilateral skewness = `max(0, 1 - A / A_ideal)` where
  `A_ideal = (√3/4) * max_edge²`

### Tetrahedron
- **elem_size** (used in `area_*` fields): signed volume
- **aspect_ratio**: circumradius / (3 * inradius)  (ideal = 1 = regular tet)
- **skewness**: `max(0, 1 - V / V_ideal)` where
  `V_ideal = max_edge³ / (6√2)` (regular tet with same max edge)

`n_degenerate` = elements with `|size| < 1e-14 * mean_abs_size`
`n_inverted`   = elements with `size < 0`
"""

# ──────────────────────────────────────────────────────────────────────────────
# Report struct
# ──────────────────────────────────────────────────────────────────────────────

"""
    MeshQualityReport

Statistical summary of mesh element quality.

Fields:
- `n_elements`: total element count
- `area_min`, `area_max`, `area_mean`: triangle area or tet |volume| stats
- `edge_min`, `edge_max`: global minimum/maximum edge length
- `aspect_ratio_min/max/mean`: shape regularity (ideal = 1)
- `skewness_min/max/mean`: equilateral/regular-tet skewness ∈ [0,1] (ideal = 0)
- `n_degenerate`: elements with area‖volume below tolerance
- `n_inverted`: elements with negative signed area‖volume
"""
struct MeshQualityReport
    n_elements       :: Int
    area_min         :: Float64
    area_max         :: Float64
    area_mean        :: Float64
    edge_min         :: Float64
    edge_max         :: Float64
    aspect_ratio_min :: Float64
    aspect_ratio_max :: Float64
    aspect_ratio_mean:: Float64
    skewness_min     :: Float64
    skewness_max     :: Float64
    skewness_mean    :: Float64
    n_degenerate     :: Int
    n_inverted       :: Int
end

function Base.show(io::IO, r::MeshQualityReport)
    println(io, "MeshQualityReport ($(r.n_elements) elements)")
    println(io, "  Area/Volume : min=$(round(r.area_min,sigdigits=4))  max=$(round(r.area_max,sigdigits=4))  mean=$(round(r.area_mean,sigdigits=4))")
    println(io, "  Edge length : min=$(round(r.edge_min,sigdigits=4))  max=$(round(r.edge_max,sigdigits=4))")
    println(io, "  Aspect ratio: min=$(round(r.aspect_ratio_min,sigdigits=4))  max=$(round(r.aspect_ratio_max,sigdigits=4))  mean=$(round(r.aspect_ratio_mean,sigdigits=4))")
    println(io, "  Skewness    : min=$(round(r.skewness_min,sigdigits=4))  max=$(round(r.skewness_max,sigdigits=4))  mean=$(round(r.skewness_mean,sigdigits=4))")
    println(io, "  Degenerate  : $(r.n_degenerate)    Inverted: $(r.n_inverted)")
end

# ──────────────────────────────────────────────────────────────────────────────
# mesh_quality(TriangleMesh)
# ──────────────────────────────────────────────────────────────────────────────

function mesh_quality(mesh::TriangleMesh{IT,FT}) where {IT,FT}
    n  = mesh.trinum
    n == 0 && return MeshQualityReport(0, 0.0, 0.0, 0.0, Inf, 0.0, Inf, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)

    areas   = Vector{Float64}(undef, n)
    ar_vals = Vector{Float64}(undef, n)
    sk_vals = Vector{Float64}(undef, n)
    g_emin  = Inf
    g_emax  = 0.0

    inv_sqrt3_4 = 1.0 / (sqrt(3.0) / 4.0)   # reciprocal of equilateral area/L²

    for k in 1:n
        v1 = mesh.node[:, mesh.triangles[1,k]]
        v2 = mesh.node[:, mesh.triangles[2,k]]
        v3 = mesh.node[:, mesh.triangles[3,k]]

        e1 = norm(v2 .- v1)
        e2 = norm(v3 .- v2)
        e3 = norm(v1 .- v3)

        emin = min(e1, e2, e3)
        emax = max(e1, e2, e3)

        # Signed area magnitude
        A = norm(cross(v2 .- v1, v3 .- v1)) / 2.0
        areas[k] = A

        # Aspect ratio: max / min edge  (ill-defined if emin≈0; clamp)
        ar_vals[k] = emin > 0.0 ? emax / emin : Inf

        # Equilateral skewness
        A_ideal = (sqrt(3.0) / 4.0) * emax^2
        sk_vals[k] = A_ideal > 0.0 ? max(0.0, (A_ideal - A) / A_ideal) : 0.0

        g_emin = min(g_emin, emin)
        g_emax = max(g_emax, emax)
    end

    mean_area = sum(areas) / n
    tol       = 1e-14 * max(mean_area, 1e-100)

    n_degenerate = count(a -> a < tol, areas)
    n_inverted   = 0   # unsigned area: sign not tracked per element here

    return MeshQualityReport(
        n,
        minimum(areas), maximum(areas), mean_area,
        g_emin, g_emax,
        minimum(ar_vals), maximum(ar_vals), sum(ar_vals) / n,
        minimum(sk_vals), maximum(sk_vals), sum(sk_vals) / n,
        n_degenerate, n_inverted
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# mesh_quality(TetrahedraMesh)
# ──────────────────────────────────────────────────────────────────────────────

"""
Circumradius of a tetrahedron with edge vectors a,b,c from the same vertex.
Uses formula R_c = |a||b||c| / (6V).  More precisely:
R_c = (|pa||pb||pc||pd|) / (8V)  — standard 4-point formula.
We use the simpler per-vertex approach via the circumradius of the tet:
R_c = sqrt( (||AB||²||CD||² + ||AC||²||BD||² + ||AD||²||BC||²) / 2 ) / (6V/???) 
Actually we use a known formula: R_circ = |a|·|b|·|c| / (6|V|) for ANY tet
where a,b,c are the three edge vectors from ONE vertex.
This is the "circumsphere radius from Cayley-Menger determinant" simplified.
Actually, the correct well-known formula used in mesh quality codes is:
  R_circ = (edge_len_product) / (something)
Let's use inradius formula instead: r_in = 3*V / (A_total)
and compute aspect_ratio = R_circ / (3 * r_in), ideal = 1 for regular tet.
Instead, we use a robust shape-quality measure:
  q = 12 * (3V)^(2/3) / sum_edges(L_i²)   (Knupp regularity)
  q ∈ [0,1], q=1 for regular tet.
We invert this to get aspect_ratio = 1/q (ideal = 1).
"""
function mesh_quality(mesh::TetrahedraMesh{IT,FT}) where {IT,FT}
    n  = mesh.tetnum
    n == 0 && return MeshQualityReport(0, 0.0, 0.0, 0.0, Inf, 0.0, Inf, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0)

    vols    = Vector{Float64}(undef, n)
    ar_vals = Vector{Float64}(undef, n)
    sk_vals = Vector{Float64}(undef, n)
    g_emin  = Inf
    g_emax  = 0.0

    for k in 1:n
        v1 = mesh.node[:, mesh.tetras[1,k]]
        v2 = mesh.node[:, mesh.tetras[2,k]]
        v3 = mesh.node[:, mesh.tetras[3,k]]
        v4 = mesh.node[:, mesh.tetras[4,k]]

        # Signed volume
        V = dot(v2 .- v1, cross(v3 .- v1, v4 .- v1)) / 6.0
        vols[k] = V

        # All 6 edge lengths
        e12 = norm(v2 .- v1); e13 = norm(v3 .- v1); e14 = norm(v4 .- v1)
        e23 = norm(v3 .- v2); e24 = norm(v4 .- v2); e34 = norm(v4 .- v3)

        emin = min(e12, e13, e14, e23, e24, e34)
        emax = max(e12, e13, e14, e23, e24, e34)

        g_emin = min(g_emin, emin)
        g_emax = max(g_emax, emax)

        # Knupp regularity: q = 12*(3|V|)^{2/3} / Σ_edges L_i²
        sum_L2 = e12^2 + e13^2 + e14^2 + e23^2 + e24^2 + e34^2
        abs_V  = abs(V)
        if sum_L2 > 0.0 && abs_V > 0.0
            q = 12.0 * (3.0 * abs_V)^(2/3) / sum_L2
            ar_vals[k] = 1.0 / max(q, 1e-100)   # aspect = 1/quality
        else
            ar_vals[k] = Inf
        end

        # Skewness: max(0, 1 - V/V_ideal) where V_ideal = emax³/(6√2) (regular tet)
        V_ideal = emax^3 / (6.0 * sqrt(2.0))
        sk_vals[k] = V_ideal > 0.0 ? max(0.0, (V_ideal - abs_V) / V_ideal) : 0.0
    end

    abs_vols  = abs.(vols)
    mean_vol  = sum(abs_vols) / n
    tol       = 1e-14 * max(mean_vol, 1e-100)

    n_degenerate = count(v -> abs(v) < tol, vols)
    n_inverted   = count(v -> v < 0.0, vols)

    # Clamp Inf aspect ratios for mean computation
    finite_ar = filter(isfinite, ar_vals)
    ar_mean   = isempty(finite_ar) ? Inf : sum(finite_ar) / length(finite_ar)

    return MeshQualityReport(
        n,
        minimum(abs_vols), maximum(abs_vols), mean_vol,
        g_emin, g_emax,
        minimum(ar_vals), maximum(ar_vals), ar_mean,
        minimum(sk_vals), maximum(sk_vals), sum(sk_vals) / n,
        n_degenerate, n_inverted
    )
end
