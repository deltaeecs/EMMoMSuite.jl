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
