using ..CoreModule

"""
    read_msh_mesh(pathname::String; FT=Float64)

Read a Gmsh (.msh) mesh file (version 4.1).
Currently supports only Triangle elements (type 2).
"""
function read_msh_mesh(pathname::String; FT = Float64)
    if !endswith(pathname, ".msh")
        error("Only .msh files are supported.")
    end

    lines = readlines(pathname)

    # Check version
    if lines[2] != "4.1 0 8"
        # Warning: strict check, might need relaxation for other 4.x versions
        # But for now let's assume 4.1 ASCII
    end

    node_start = findfirst(==("\$Nodes"), lines)
    elem_start = findfirst(==("\$Elements"), lines)

    if isnothing(node_start) || isnothing(elem_start)
        error("Invalid Gmsh file structure.")
    end

    # Parse Nodes
    # $Nodes
    # numEntityBlocks(1) numNodes(1) minNodeTag(1) maxNodeTag(1)
    # entityDim(1) entityTag(1) parametric(1) numNodesInBlock(1)
    # nodeTag(1)
    # ...
    # x(1) y(1) z(1)
    # ...
    # $EndNodes

    # Simplified parser: just look for blocks and read nodes
    # We need to map nodeTag to index 1..N

    # Read header of Nodes section
    parts = split(lines[node_start+1])
    num_blocks = parse(Int, parts[1])
    total_nodes = parse(Int, parts[2])

    node = zeros(FT, 3, total_nodes)
    node_tag_map = Dict{Int,Int}()

    current_line = node_start + 2
    global_node_idx = 1

    for b = 1:num_blocks
        parts = split(lines[current_line])
        # entityDim = parse(Int, parts[1])
        # entityTag = parse(Int, parts[2])
        # parametric = parse(Int, parts[3])
        num_nodes_in_block = parse(Int, parts[4])
        current_line += 1

        # Read tags
        tags = Int[]
        for i = 1:num_nodes_in_block
            push!(tags, parse(Int, lines[current_line]))
            current_line += 1
        end

        # Read coordinates
        for i = 1:num_nodes_in_block
            coords = parse.(FT, split(lines[current_line]))
            node[:, global_node_idx] = coords[1:3] # x, y, z
            node_tag_map[tags[i]] = global_node_idx
            global_node_idx += 1
            current_line += 1
        end
    end

    # Parse Elements
    # $Elements
    # numEntityBlocks(1) numElements(1) minElementTag(1) maxElementTag(1)
    # entityDim(1) entityTag(1) elementType(1) numElementsInBlock(1)
    # elementTag(1) nodeTag(1) ... nodeTag(k)
    # ...
    # $EndElements

    parts = split(lines[elem_start+1])
    num_elem_blocks = parse(Int, parts[1])
    total_elems = parse(Int, parts[2])

    # We only care about triangles (type 2)
    # First pass to count triangles
    num_triangles = 0
    temp_line = elem_start + 2

    block_infos = [] # Store (line_idx, num_elems, type, entityTag)

    for b = 1:num_elem_blocks
        parts = split(lines[temp_line])
        entity_dim = parse(Int, parts[1])
        entity_tag = parse(Int, parts[2])
        elem_type = parse(Int, parts[3])
        num_elems_in_block = parse(Int, parts[4])

        push!(block_infos, (temp_line + 1, num_elems_in_block, elem_type, entity_tag))

        if elem_type == 2 # 3-node triangle
            num_triangles += num_elems_in_block
        end

        temp_line += 1 + num_elems_in_block
    end

    triangles = zeros(Int, 3, num_triangles)
    tags = zeros(Int, num_triangles)

    tri_idx = 1
    for (start_line, num, type, tag) in block_infos
        if type == 2
            for i = 0:num-1
                parts = split(lines[start_line+i])
                # elementTag node1 node2 node3
                # n1 = parse(Int, parts[2])
                # n2 = parse(Int, parts[3])
                # n3 = parse(Int, parts[4])

                # Gmsh 4.1: elementTag nodeTags...
                n1 = node_tag_map[parse(Int, parts[2])]
                n2 = node_tag_map[parse(Int, parts[3])]
                n3 = node_tag_map[parse(Int, parts[4])]

                triangles[:, tri_idx] = [n1, n2, n3]
                tags[tri_idx] = tag # Use entityTag as physical tag/region ID
                tri_idx += 1
            end
        end
    end

    return TriangleMesh(num_triangles, node, triangles, tags)
end
