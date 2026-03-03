using StaticArrays
using LinearAlgebra

"""
    generate_rectangle_mesh(Lx, Ly, nx, ny)

Generate a rectangular plate mesh in the xy-plane.
Lx, Ly: Dimensions
nx, ny: Number of segments along x and y
"""
function generate_rectangle_mesh(Lx::FT, Ly::FT, nx::Int, ny::Int) where {FT<:Real}
    dx = Lx / nx
    dy = Ly / ny

    # Nodes
    num_nodes = (nx + 1) * (ny + 1)
    nodes = zeros(FT, 3, num_nodes)

    idx = 1
    for j = 0:ny
        y = j * dy - Ly / 2
        for i = 0:nx
            x = i * dx - Lx / 2
            nodes[1, idx] = x
            nodes[2, idx] = y
            nodes[3, idx] = 0.0
            idx += 1
        end
    end

    # Elements (Triangles)
    # Each rectangle (i,j) is split into 2 triangles
    num_elems = 2 * nx * ny
    elements = zeros(Int, 3, num_elems)

    elem_idx = 1
    for j = 1:ny
        for i = 1:nx
            # Node indices (1-based)
            n1 = (j - 1) * (nx + 1) + i
            n2 = n1 + 1
            n3 = j * (nx + 1) + i
            n4 = n3 + 1

            # Triangle 1 (n1, n2, n3)
            elements[1, elem_idx] = n1
            elements[2, elem_idx] = n2
            elements[3, elem_idx] = n3
            elem_idx += 1

            # Triangle 2 (n2, n4, n3)
            elements[1, elem_idx] = n2
            elements[2, elem_idx] = n4
            elements[3, elem_idx] = n3
            elem_idx += 1
        end
    end

    return TriangleMesh(num_elems, nodes, elements)
end

"""
    generate_cylinder_mesh(radius, height, n_circum, n_height; closed=true)

Generate a cylinder mesh aligned with z-axis.
"""
function generate_cylinder_mesh(
    radius::FT,
    height::FT,
    n_circum::Int,
    n_height::Int;
    closed::Bool = true,
) where {FT<:Real}
    # Side nodes
    # (n_circum) points per circle (last = first) -> n_circum unique points
    # (n_height + 1) layers

    num_side_nodes = n_circum * (n_height + 1)
    nodes = zeros(FT, 3, num_side_nodes) # Will resize later if closed

    dphi = 2π / n_circum
    dz = height / n_height

    idx = 1
    for j = 0:n_height
        z = j * dz - height / 2
        for i = 1:n_circum
            phi = (i - 1) * dphi
            nodes[1, idx] = radius * cos(phi)
            nodes[2, idx] = radius * sin(phi)
            nodes[3, idx] = z
            idx += 1
        end
    end

    # Side elements
    num_side_elems = 2 * n_circum * n_height
    elements = zeros(Int, 3, num_side_elems) # Will resize later

    elem_idx = 1
    for j = 1:n_height
        for i = 1:n_circum
            # Current layer
            n1 = (j - 1) * n_circum + i
            n2 = (i == n_circum) ? (j - 1) * n_circum + 1 : n1 + 1

            # Next layer
            n3 = n1 + n_circum
            n4 = n2 + n_circum

            # Tri 1
            elements[1, elem_idx] = n1
            elements[2, elem_idx] = n2
            elements[3, elem_idx] = n3
            elem_idx += 1

            # Tri 2
            elements[1, elem_idx] = n2
            elements[2, elem_idx] = n4
            elements[3, elem_idx] = n3
            elem_idx += 1
        end
    end

    if closed
        # Add top and bottom caps
        # Simple fan triangulation with center point

        # Add center nodes
        nodes = hcat(nodes, [0.0, 0.0, -height / 2]) # Bottom center
        bottom_center_idx = size(nodes, 2)

        nodes = hcat(nodes, [0.0, 0.0, height / 2]) # Top center
        top_center_idx = size(nodes, 2)

        # Bottom cap
        bottom_elems = zeros(Int, 3, n_circum)
        for i = 1:n_circum
            n1 = i
            n2 = (i == n_circum) ? 1 : i + 1
            # Clockwise looking from outside (bottom normal points -z)
            # Normal should be outward. 
            # Side normals are outward.
            # Bottom normal is -z.
            # Triangle (n1, n2, center) -> (x1, y1) x (x2, y2) ~ z
            # We want -z. So (n2, n1, center)
            bottom_elems[1, i] = n2
            bottom_elems[2, i] = n1
            bottom_elems[3, i] = bottom_center_idx
        end

        # Top cap
        top_elems = zeros(Int, 3, n_circum)
        start_idx = n_circum * n_height
        for i = 1:n_circum
            n1 = start_idx + i
            n2 = (i == n_circum) ? start_idx + 1 : n1 + 1
            # Top normal is +z. (n1, n2, center)
            top_elems[1, i] = n1
            top_elems[2, i] = n2
            top_elems[3, i] = top_center_idx
        end

        elements = hcat(elements, bottom_elems, top_elems)
    end

    return TriangleMesh(size(elements, 2), nodes, elements)
end

"""
    generate_sphere_mesh(radius, n_theta, n_phi)

Generate a UV sphere mesh.
radius: Sphere radius
n_theta: Number of segments along theta (latitude).
n_phi: Number of segments along phi (longitude).
"""
function generate_sphere_mesh(radius::FT, n_theta::Int, n_phi::Int) where {FT<:Real}
    # Vertices
    # North Pole + (n_theta - 1) rings * n_phi + South Pole
    num_nodes = 1 + (n_theta - 1) * n_phi + 1
    nodes = zeros(FT, 3, num_nodes)

    # North Pole
    nodes[1, 1] = 0.0
    nodes[2, 1] = 0.0
    nodes[3, 1] = radius

    # Rings
    idx = 2
    for i = 1:(n_theta-1)
        theta = i * FT(π) / n_theta
        sin_theta = sin(theta)
        cos_theta = cos(theta)

        for j = 0:(n_phi-1)
            phi = j * 2 * FT(π) / n_phi
            nodes[1, idx] = radius * sin_theta * cos(phi)
            nodes[2, idx] = radius * sin_theta * sin(phi)
            nodes[3, idx] = radius * cos_theta
            idx += 1
        end
    end

    # South Pole
    nodes[1, num_nodes] = 0.0
    nodes[2, num_nodes] = 0.0
    nodes[3, num_nodes] = -radius

    # Elements
    # Top cap: n_phi
    # Bands: (n_theta - 2) * 2 * n_phi
    # Bottom cap: n_phi
    num_elems = 2 * n_phi * (n_theta - 1)
    elements = zeros(Int, 3, num_elems)

    elem_idx = 1

    # Top Cap
    # North Pole is 1
    # First ring starts at 2
    for j = 0:(n_phi-1)
        n1 = 1
        n2 = 2 + j
        n3 = 2 + (j + 1) % n_phi

        elements[1, elem_idx] = n1
        elements[2, elem_idx] = n2
        elements[3, elem_idx] = n3
        elem_idx += 1
    end

    # Middle Bands
    for i = 1:(n_theta-2)
        s_curr = 2 + (i - 1) * n_phi
        s_next = s_curr + n_phi

        for j = 0:(n_phi-1)
            p1 = s_curr + j
            p2 = s_curr + (j + 1) % n_phi
            p3 = s_next + j
            p4 = s_next + (j + 1) % n_phi

            # Triangle 1
            elements[1, elem_idx] = p1
            elements[2, elem_idx] = p3
            elements[3, elem_idx] = p2
            elem_idx += 1

            # Triangle 2
            elements[1, elem_idx] = p2
            elements[2, elem_idx] = p3
            elements[3, elem_idx] = p4
            elem_idx += 1
        end
    end

    # Bottom Cap
    s_last = 2 + (n_theta - 2) * n_phi
    n_sp = num_nodes

    for j = 0:(n_phi-1)
        p1 = s_last + j
        p2 = s_last + (j + 1) % n_phi
        p3 = n_sp

        # Triangle (p1, p3, p2) for outward normal
        elements[1, elem_idx] = p1
        elements[2, elem_idx] = p3
        elements[3, elem_idx] = p2
        elem_idx += 1
    end

    return TriangleMesh(num_elems, nodes, elements)
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 16.1 — 新增表面网格生成器
# ─────────────────────────────────────────────────────────────────────────────

"""
    generate_ellipsoid_mesh(a, b, c, n_theta, n_phi)

Generate an ellipsoidal surface mesh with semi-axes a (x), b (y), c (z).
Uses the same UV topology as `generate_sphere_mesh`.
"""
function generate_ellipsoid_mesh(
    a::FT, b::FT, c::FT, n_theta::Int, n_phi::Int
) where {FT<:Real}
    sphere = generate_sphere_mesh(one(FT), n_theta, n_phi)
    nodes = copy(sphere.node)
    nodes[1, :] .*= a
    nodes[2, :] .*= b
    nodes[3, :] .*= c
    return TriangleMesh(sphere.trinum, nodes, sphere.triangles, sphere.tags)
end

"""
    generate_cone_mesh(radius, height, n_circum, n_height; closed=true)

Generate a cone mesh with apex at +z/2 and circular base at -z/2.
`n_height` is the number of lateral ring-bands (≥1).
When `closed=true`, a base cap is added.
"""
function generate_cone_mesh(
    radius::FT,
    height::FT,
    n_circum::Int,
    n_height::Int;
    closed::Bool = true,
) where {FT<:Real}
    n_rings = n_height            # rings 0..n_height-1 with n_circum nodes each
    num_ring_nodes = n_rings * n_circum
    # apex is one extra node
    all_nodes = zeros(FT, 3, num_ring_nodes + 1 + (closed ? 1 : 0))
    for i in 0:(n_rings-1)
        t  = FT(i) / FT(n_height)
        r_i = radius * (one(FT) - t)
        z_i = -height / 2 + t * height
        for j in 0:(n_circum-1)
            phi = j * 2 * FT(π) / n_circum
            idx = i * n_circum + j + 1
            all_nodes[1, idx] = r_i * cos(phi)
            all_nodes[2, idx] = r_i * sin(phi)
            all_nodes[3, idx] = z_i
        end
    end
    apex_idx = num_ring_nodes + 1
    all_nodes[3, apex_idx] = height / 2

    n_lateral_full = max(0, (n_rings - 1)) * 2 * n_circum
    n_apex_band    = n_circum
    n_cap          = closed ? n_circum : 0
    num_elems      = n_lateral_full + n_apex_band + n_cap

    elements = zeros(Int, 3, num_elems)
    eidx = 1

    # Full bands between consecutive rings
    for i in 0:(n_rings-2)
        s_curr = i * n_circum + 1
        s_next = (i + 1) * n_circum + 1
        for j in 0:(n_circum-1)
            p1 = s_curr + j;           p2 = s_curr + (j + 1) % n_circum
            p3 = s_next + j;           p4 = s_next + (j + 1) % n_circum
            elements[:, eidx] = [p1, p3, p2]; eidx += 1
            elements[:, eidx] = [p2, p3, p4]; eidx += 1
        end
    end

    # Apex band (last ring → apex)
    s_last = (n_rings - 1) * n_circum + 1
    for j in 0:(n_circum-1)
        p1 = s_last + j
        p2 = s_last + (j + 1) % n_circum
        elements[:, eidx] = [p1, apex_idx, p2]; eidx += 1
    end

    # Base cap, outward normal −z: CCW when viewed from −z side
    if closed
        base_idx = apex_idx + 1
        all_nodes[3, base_idx] = -height / 2   # already pre-allocated
        for j in 0:(n_circum-1)
            p1 = 1 + j
            p2 = 1 + (j + 1) % n_circum
            elements[:, eidx] = [p2, p1, base_idx]; eidx += 1
        end
    end

    all_nodes = all_nodes[:, 1:(apex_idx + (closed ? 1 : 0))]
    return TriangleMesh(num_elems, all_nodes, elements)
end

"""
    generate_torus_mesh(R, r, n_major, n_minor)

Generate a torus surface mesh.
`R` = major radius (centre-line), `r` = tube radius.
`n_major` segments around the z-axis, `n_minor` around the tube.
"""
function generate_torus_mesh(
    R::FT, r::FT, n_major::Int, n_minor::Int
) where {FT<:Real}
    num_nodes = n_major * n_minor
    nodes = zeros(FT, 3, num_nodes)
    for i in 0:(n_major-1)
        theta = i * 2 * FT(π) / n_major
        ct, st = cos(theta), sin(theta)
        for j in 0:(n_minor-1)
            phi = j * 2 * FT(π) / n_minor
            cp, sp = cos(phi), sin(phi)
            idx = i * n_minor + j + 1
            nodes[1, idx] = (R + r * cp) * ct
            nodes[2, idx] = (R + r * cp) * st
            nodes[3, idx] = r * sp
        end
    end

    num_elems = 2 * n_major * n_minor
    elements = zeros(Int, 3, num_elems)
    eidx = 1
    for i in 0:(n_major-1)
        ni      = i * n_minor + 1
        ni_next = ((i + 1) % n_major) * n_minor + 1
        for j in 0:(n_minor-1)
            p1 = ni + j;           p2 = ni + (j + 1) % n_minor
            p3 = ni_next + j;      p4 = ni_next + (j + 1) % n_minor
            elements[:, eidx] = [p1, p3, p2]; eidx += 1
            elements[:, eidx] = [p2, p3, p4]; eidx += 1
        end
    end

    return TriangleMesh(num_elems, nodes, elements)
end

# ─────────────────────────────────────────────────────────────────────────────
# Phase 16.1 — 体网格生成器
# ─────────────────────────────────────────────────────────────────────────────

"""
    generate_box_volume_mesh(Lx, Ly, Lz, nx, ny, nz)

Generate a structured hexahedral body mesh for a box centred at the origin.
Returns a `HexahedraMesh` with `nx*ny*nz` elements.
Node ordering per hex (Nastran CHEXA convention):
  nodes 1–4 = bottom face (k), CCW viewed from −z;
  nodes 5–8 = top face (k+1), same plan order.
"""
function generate_box_volume_mesh(
    Lx::FT, Ly::FT, Lz::FT, nx::Int, ny::Int, nz::Int
) where {FT<:Real}
    dx, dy, dz = Lx / nx, Ly / ny, Lz / nz

    nidx(i, j, k) = i + (nx + 1) * j + (nx + 1) * (ny + 1) * k + 1
    num_nodes = (nx + 1) * (ny + 1) * (nz + 1)

    nodes = zeros(FT, 3, num_nodes)
    for k in 0:nz, j in 0:ny, i in 0:nx
        idx = nidx(i, j, k)
        nodes[1, idx] = i * dx - Lx / 2
        nodes[2, idx] = j * dy - Ly / 2
        nodes[3, idx] = k * dz - Lz / 2
    end

    num_elems = nx * ny * nz
    hexes = zeros(Int, 8, num_elems)
    eidx = 1
    for k in 0:(nz-1), j in 0:(ny-1), i in 0:(nx-1)
        hexes[1, eidx] = nidx(i,   j,   k)
        hexes[2, eidx] = nidx(i+1, j,   k)
        hexes[3, eidx] = nidx(i+1, j+1, k)
        hexes[4, eidx] = nidx(i,   j+1, k)
        hexes[5, eidx] = nidx(i,   j,   k+1)
        hexes[6, eidx] = nidx(i+1, j,   k+1)
        hexes[7, eidx] = nidx(i+1, j+1, k+1)
        hexes[8, eidx] = nidx(i,   j+1, k+1)
        eidx += 1
    end

    return HexahedraMesh(num_elems, nodes, hexes)
end

"""
    generate_box_tet_mesh(Lx, Ly, Lz, nx, ny, nz)

Generate a structured tetrahedral body mesh for a box centred at the origin.
Each hexahedral cell is split into 6 tetrahedra using the consistent
v[1]–v[7] main diagonal, giving `6*nx*ny*nz` elements.
"""
function generate_box_tet_mesh(
    Lx::FT, Ly::FT, Lz::FT, nx::Int, ny::Int, nz::Int
) where {FT<:Real}
    dx, dy, dz = Lx / nx, Ly / ny, Lz / nz

    nidx(i, j, k) = i + (nx + 1) * j + (nx + 1) * (ny + 1) * k + 1
    num_nodes = (nx + 1) * (ny + 1) * (nz + 1)

    nodes = zeros(FT, 3, num_nodes)
    for k in 0:nz, j in 0:ny, i in 0:nx
        idx = nidx(i, j, k)
        nodes[1, idx] = i * dx - Lx / 2
        nodes[2, idx] = j * dy - Ly / 2
        nodes[3, idx] = k * dz - Lz / 2
    end

    num_elems = 6 * nx * ny * nz
    tetras = zeros(Int, 4, num_elems)
    eidx = 1
    for k in 0:(nz-1), j in 0:(ny-1), i in 0:(nx-1)
        v = (
            nidx(i,   j,   k),    # v1 = (0,0,0)
            nidx(i+1, j,   k),    # v2 = (1,0,0)
            nidx(i+1, j+1, k),    # v3 = (1,1,0)
            nidx(i,   j+1, k),    # v4 = (0,1,0)
            nidx(i,   j,   k+1),  # v5 = (0,0,1)
            nidx(i+1, j,   k+1),  # v6 = (1,0,1)
            nidx(i+1, j+1, k+1),  # v7 = (1,1,1)
            nidx(i,   j+1, k+1),  # v8 = (0,1,1)
        )
        # Freudenthal triangulation (6 tets, all with positive signed volume)
        # diagonal v1→v7; each tet verified: tet_volume > 0
        tetras[:, eidx] = [v[1], v[2], v[3], v[7]]; eidx += 1  # σ=(1,2,3)
        tetras[:, eidx] = [v[1], v[2], v[7], v[6]]; eidx += 1  # σ=(1,3,2)*
        tetras[:, eidx] = [v[1], v[4], v[7], v[3]]; eidx += 1  # σ=(2,1,3)*
        tetras[:, eidx] = [v[1], v[4], v[8], v[7]]; eidx += 1  # σ=(2,3,1)
        tetras[:, eidx] = [v[1], v[5], v[6], v[7]]; eidx += 1  # σ=(3,1,2)
        tetras[:, eidx] = [v[1], v[5], v[7], v[8]]; eidx += 1  # σ=(3,2,1)*
    end

    return TetrahedraMesh(num_elems, nodes, tetras)
end
