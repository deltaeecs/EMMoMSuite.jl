using ..CoreModule

function parse_nastran_field(line, start_col, end_col)
    if length(line) < start_col
        return ""
    end
    e = min(length(line), end_col)
    return strip(line[start_col:e])
end

"""
    parse_nastran_float(s::AbstractString, ::Type{FT}) where FT

Parse a Nastran float string, handling compressed scientific notation.
Nastran allows formats like `1.78-15` meaning `1.78e-15`, and `3.2+5` meaning `3.2e+5`.
"""
function parse_nastran_float(s::AbstractString, ::Type{FT}) where FT
    s = strip(s)
    isempty(s) && return zero(FT)
    # Try normal parse first
    v = tryparse(FT, s)
    !isnothing(v) && return v
    # Handle Nastran compressed scientific notation: e.g. "1.78-15" -> "1.78e-15"
    # Pattern: digits with optional dot, followed by +/- and digits (no E)
    m = match(r"^([+-]?\d*\.?\d+)([+-]\d+)$", s)
    if !isnothing(m)
        return parse(FT, m.captures[1] * "e" * m.captures[2])
    end
    # Fallback
    return parse(FT, s)
end

"""
    read_nas_mesh(pathname::String; FT=Float64, scale=1.0)

Read a Nastran (.nas) mesh file.
"""
function read_nas_mesh(pathname::String; FT=Float64, scale=1.0)
    if !endswith(pathname, ".nas")
        error("Only .nas files are supported.")
    end

    linenum, nodenum, trinum, tetranum, hexanum = open(pathname, "r") do f
        linenum = 0; nodenum = 0; trinum = 0; tetranum = 0; hexanum = 0
        while !eof(f)
            linenum += 1
            line = readline(f)
            startswith(line, r"\$|//|#") && continue
            if occursin("GRID", line)
                nodenum += 1
            elseif occursin("CTRIA3", line)
                trinum += 1
            elseif occursin("CTETRA", line)
                tetranum += 1
            elseif occursin("CHEXA", line)
                hexanum += 1
            end
        end
        linenum, nodenum, trinum, tetranum, hexanum
    end

    node = zeros(FT, 3, nodenum)
    triangles = zeros(Int, 3, trinum)
    tetras = zeros(Int, 4, tetranum)
    hexas = zeros(Int, 8, hexanum)
    tri_tags = zeros(Int, trinum)
    tet_tags = zeros(Int, tetranum)
    hex_tags = zeros(Int, hexanum)
    
    # Second pass: read data
    tri_idx, tet_idx, hex_idx = open(pathname, "r") do f
        node_idx = 1
        tri_idx = 1
        tet_idx = 1
        hex_idx = 1
        
        # Map from NAS grid ID to 1-based index
        grid_id_map = Dict{Int, Int}()
        
        while !eof(f)
            line = readline(f)
            startswith(line, r"\$|//|#") && continue
            
            if startswith(line, "GRID")
                if occursin(",", line)
                    parts = split(line, ",")
                    parts = strip.(parts)
                    if length(parts) >= 6
                        id = parse(Int, parts[2])
                        x = parse_nastran_float(parts[4], FT)
                        y = parse_nastran_float(parts[5], FT)
                        z = parse_nastran_float(parts[6], FT)
                        
                        node[:, node_idx] = [x, y, z] .* scale
                        grid_id_map[id] = node_idx
                        node_idx += 1
                    end
                else
                    # Fixed width (Small Field)
                    if startswith(line, "GRID*")
                        # Large Field Format
                        # Line 1: GRID* (1-8), ID (9-24), CP (25-40), X (41-56), Y (57-72)
                        id_str = parse_nastran_field(line, 9, 24)
                        x_str = parse_nastran_field(line, 41, 56)
                        y_str = parse_nastran_field(line, 57, 72)
                        
                        # Read next line for Z
                        if !eof(f)
                            line2 = readline(f)
                            # Line 2: * (1-8), Z (9-24)
                            z_str = parse_nastran_field(line2, 9, 24)
                            
                            if !isempty(id_str) && !isempty(x_str) && !isempty(y_str) && !isempty(z_str)
                                id = parse(Int, id_str)
                                x = parse_nastran_float(x_str, FT)
                                y = parse_nastran_float(y_str, FT)
                                z = parse_nastran_float(z_str, FT)
                                
                                node[:, node_idx] = [x, y, z] .* scale
                                grid_id_map[id] = node_idx
                                node_idx += 1
                            end
                        end
                        continue
                    end
                    
                    id_str = parse_nastran_field(line, 9, 16)
                    x_str = parse_nastran_field(line, 25, 32)
                    y_str = parse_nastran_field(line, 33, 40)
                    z_str = parse_nastran_field(line, 41, 48)
                    
                    if !isempty(id_str) && !isempty(x_str) && !isempty(y_str) && !isempty(z_str)
                        id = parse(Int, id_str)
                        x = parse_nastran_float(x_str, FT)
                        y = parse_nastran_float(y_str, FT)
                        z = parse_nastran_float(z_str, FT)
                        
                        node[:, node_idx] = [x, y, z] .* scale
                        grid_id_map[id] = node_idx
                        node_idx += 1
                    end
                end
            elseif startswith(line, "CTRIA3")
                if occursin(",", line)
                    parts = split(line, ",")
                    parts = strip.(parts)
                    if length(parts) >= 6
                        pid = parse(Int, parts[3])
                        g1 = parse(Int, parts[4])
                        g2 = parse(Int, parts[5])
                        g3 = parse(Int, parts[6])
                        
                        if haskey(grid_id_map, g1) && haskey(grid_id_map, g2) && haskey(grid_id_map, g3)
                            triangles[:, tri_idx] = [grid_id_map[g1], grid_id_map[g2], grid_id_map[g3]]
                            tri_tags[tri_idx] = pid
                            tri_idx += 1
                        end
                    end
                else
                    # Fixed width
                    pid_str = parse_nastran_field(line, 17, 24)
                    g1_str = parse_nastran_field(line, 25, 32)
                    g2_str = parse_nastran_field(line, 33, 40)
                    g3_str = parse_nastran_field(line, 41, 48)
                    
                    if !isempty(pid_str) && !isempty(g1_str) && !isempty(g2_str) && !isempty(g3_str)
                        pid = parse(Int, pid_str)
                        g1 = parse(Int, g1_str)
                        g2 = parse(Int, g2_str)
                        g3 = parse(Int, g3_str)
                        
                        if haskey(grid_id_map, g1) && haskey(grid_id_map, g2) && haskey(grid_id_map, g3)
                            triangles[:, tri_idx] = [grid_id_map[g1], grid_id_map[g2], grid_id_map[g3]]
                            tri_tags[tri_idx] = pid
                            tri_idx += 1
                        end
                    end
                end
            elseif startswith(line, "CTETRA")
                if occursin(",", line)
                    parts = split(line, ",")
                    parts = strip.(parts)
                    if length(parts) >= 7
                        # CTETRA, EID, PID, G1, G2, G3, G4
                        pid = parse(Int, parts[3])
                        g1 = parse(Int, parts[4])
                        g2 = parse(Int, parts[5])
                        g3 = parse(Int, parts[6])
                        g4 = parse(Int, parts[7])
                        
                        if haskey(grid_id_map, g1) && haskey(grid_id_map, g2) && haskey(grid_id_map, g3) && haskey(grid_id_map, g4)
                            tetras[:, tet_idx] = [grid_id_map[g1], grid_id_map[g2], grid_id_map[g3], grid_id_map[g4]]
                            tet_tags[tet_idx] = pid
                            tet_idx += 1
                        end
                    end
                else
                    # Fixed width
                    pid_str = parse_nastran_field(line, 17, 24)
                    g1_str = parse_nastran_field(line, 25, 32)
                    g2_str = parse_nastran_field(line, 33, 40)
                    g3_str = parse_nastran_field(line, 41, 48)
                    g4_str = parse_nastran_field(line, 49, 56)
                    
                    if !isempty(pid_str) && !isempty(g1_str) && !isempty(g2_str) && !isempty(g3_str) && !isempty(g4_str)
                        pid = parse(Int, pid_str)
                        g1 = parse(Int, g1_str)
                        g2 = parse(Int, g2_str)
                        g3 = parse(Int, g3_str)
                        g4 = parse(Int, g4_str)
                        
                        if haskey(grid_id_map, g1) && haskey(grid_id_map, g2) && haskey(grid_id_map, g3) && haskey(grid_id_map, g4)
                            tetras[:, tet_idx] = [grid_id_map[g1], grid_id_map[g2], grid_id_map[g3], grid_id_map[g4]]
                            tet_tags[tet_idx] = pid
                            tet_idx += 1
                        end
                    end
                end
            elseif startswith(line, "CHEXA")
                hex_nodes = _parse_chexa(line, f, grid_id_map)
                if !isnothing(hex_nodes)
                    hexas[:, hex_idx] = hex_nodes[1]
                    hex_tags[hex_idx] = hex_nodes[2]
                    hex_idx += 1
                end
            end
        end
        tri_idx, tet_idx, hex_idx
    end
    
    # Resize arrays if some elements were skipped
    if tri_idx <= trinum
        triangles = triangles[:, 1:tri_idx-1]
        tri_tags = tri_tags[1:tri_idx-1]
        trinum = tri_idx - 1
    end
    
    if tet_idx <= tetranum
        tetras = tetras[:, 1:tet_idx-1]
        tet_tags = tet_tags[1:tet_idx-1]
        tetranum = tet_idx - 1
    end

    if hex_idx <= hexanum
        hexas = hexas[:, 1:hex_idx-1]
        hex_tags = hex_tags[1:hex_idx-1]
        hexanum = hex_idx - 1
    end

    if hexanum > 0 && trinum == 0 && tetranum == 0
        return HexahedraMesh(hexanum, node, hexas, hex_tags)
    elseif tetranum > 0
        return TetrahedraMesh(tetranum, node, tetras, tet_tags)
    else
        return TriangleMesh(trinum, node, triangles, tri_tags)
    end
end

"""
    _parse_chexa(line, f, grid_id_map)

Parse a CHEXA element from Nastran format (handles continuation lines).
Returns `(node_ids::Vector{Int}, pid::Int)` or `nothing` on failure.
"""
function _parse_chexa(line::String, f::IO, grid_id_map::Dict{Int,Int})
    if occursin(",", line)
        # Free-field format
        parts = split(line, ",")
        parts = strip.(parts)
        # CHEXA, EID, PID, G1..G6 (may have continuation)
        # If <=8 fields, read continuation
        all_parts = copy(parts)
        while !eof(f) && length(all_parts) < 11  # Need at least CHEXA,EID,PID,G1-G8
            next_line = readline(f)
            if startswith(next_line, "+") || startswith(next_line, ",")
                extra = split(next_line, ",")
                append!(all_parts, strip.(extra))
            else
                break
            end
        end
        pid = parse(Int, all_parts[3])
        gnodes = Int[]
        for k in 4:min(length(all_parts), 11)
            s = all_parts[k]
            isempty(s) && continue
            startswith(s, "+") && continue
            push!(gnodes, parse(Int, s))
        end
        if length(gnodes) >= 8
            mapped = [get(grid_id_map, g, 0) for g in gnodes[1:8]]
            all(>(0), mapped) && return (mapped, pid)
        end
    else
        # Fixed-width format: CHEXA has 8 fields per line
        # Line 1: CHEXA(1-8) EID(9-16) PID(17-24) G1(25-32) G2(33-40) G3(41-48) G4(49-56) G5(57-64) G6(65-72)+cont
        # Line 2: +cont(1-8) G7(9-16) G8(17-24)
        pid_str = parse_nastran_field(line, 17, 24)
        g_strs = [parse_nastran_field(line, 25 + 8*(k-1), 32 + 8*(k-1)) for k in 1:6]
        
        # Read continuation line for G7, G8
        if !eof(f)
            cont_line = readline(f)
            if startswith(cont_line, "+") || startswith(strip(cont_line), "+") || !startswith(cont_line, r"\$|GRID|CTRIA|CTETRA|CHEXA|END|BEGIN")
                push!(g_strs, parse_nastran_field(cont_line, 9, 16))
                push!(g_strs, parse_nastran_field(cont_line, 17, 24))
            end
        end
        
        if !isempty(pid_str) && length(g_strs) >= 8 && all(!isempty, g_strs[1:8])
            pid = parse(Int, pid_str)
            gnodes = [parse(Int, s) for s in g_strs[1:8]]
            mapped = [get(grid_id_map, g, 0) for g in gnodes]
            all(>(0), mapped) && return (mapped, pid)
        end
    end
    return nothing
end

"""
    read_mixed_nas_mesh(pathname::String; FT=Float64, scale=1.0)

Read a Nastran (.nas) mesh file and return surface, tetrahedral, and hexahedral meshes.
Returns `(surface_mesh::TriangleMesh, tetra_mesh::TetrahedraMesh, hexa_mesh::HexahedraMesh)`.
"""
function read_mixed_nas_mesh(pathname::String; FT=Float64, scale=1.0)
    if !endswith(pathname, ".nas")
        error("Only .nas files are supported.")
    end

    linenum, nodenum, trinum, tetranum, hexanum = open(pathname, "r") do f
        linenum = 0; nodenum = 0; trinum = 0; tetranum = 0; hexanum = 0
        while !eof(f)
            linenum += 1
            line = readline(f)
            startswith(line, r"\$|//|#") && continue
            if occursin("GRID", line)
                nodenum += 1
            elseif occursin("CTRIA3", line)
                trinum += 1
            elseif occursin("CTETRA", line)
                tetranum += 1
            elseif occursin("CHEXA", line)
                hexanum += 1
            end
        end
        linenum, nodenum, trinum, tetranum, hexanum
    end

    node = zeros(FT, 3, nodenum)
    triangles = zeros(Int, 3, trinum)
    tetras = zeros(Int, 4, tetranum)
    hexas = zeros(Int, 8, hexanum)
    tri_tags = zeros(Int, trinum)
    tet_tags = zeros(Int, tetranum)
    hex_tags = zeros(Int, hexanum)
    
    # Second pass: read data
    tri_idx, tet_idx, hex_idx = open(pathname, "r") do f
        node_idx = 1
        tri_idx = 1
        tet_idx = 1
        hex_idx = 1
        
        # Map from NAS grid ID to 1-based index
        grid_id_map = Dict{Int, Int}()
        
        while !eof(f)
            line = readline(f)
            startswith(line, r"\$|//|#") && continue
            
            if startswith(line, "GRID")
                if occursin(",", line)
                    parts = split(line, ",")
                    parts = strip.(parts)
                    if length(parts) >= 6
                        id = parse(Int, parts[2])
                        x = parse_nastran_float(parts[4], FT)
                        y = parse_nastran_float(parts[5], FT)
                        z = parse_nastran_float(parts[6], FT)
                        
                        node[:, node_idx] = [x, y, z] .* scale
                        grid_id_map[id] = node_idx
                        node_idx += 1
                    end
                else
                    # Fixed width (Small Field)
                    if startswith(line, "GRID*")
                        # Large Field Format
                        # Line 1: GRID* (1-8), ID (9-24), CP (25-40), X (41-56), Y (57-72)
                        id_str = parse_nastran_field(line, 9, 24)
                        x_str = parse_nastran_field(line, 41, 56)
                        y_str = parse_nastran_field(line, 57, 72)
                        
                        # Read next line for Z
                        if !eof(f)
                            line2 = readline(f)
                            # Line 2: * (1-8), Z (9-24)
                            z_str = parse_nastran_field(line2, 9, 24)
                            
                            if !isempty(id_str) && !isempty(x_str) && !isempty(y_str) && !isempty(z_str)
                                id = parse(Int, id_str)
                                x = parse_nastran_float(x_str, FT)
                                y = parse_nastran_float(y_str, FT)
                                z = parse_nastran_float(z_str, FT)
                                
                                node[:, node_idx] = [x, y, z] .* scale
                                grid_id_map[id] = node_idx
                                node_idx += 1
                            end
                        end
                        continue
                    end
                    
                    id_str = parse_nastran_field(line, 9, 16)
                    x_str = parse_nastran_field(line, 25, 32)
                    y_str = parse_nastran_field(line, 33, 40)
                    z_str = parse_nastran_field(line, 41, 48)
                    
                    if !isempty(id_str) && !isempty(x_str) && !isempty(y_str) && !isempty(z_str)
                        id = parse(Int, id_str)
                        x = parse_nastran_float(x_str, FT)
                        y = parse_nastran_float(y_str, FT)
                        z = parse_nastran_float(z_str, FT)
                        
                        node[:, node_idx] = [x, y, z] .* scale
                        grid_id_map[id] = node_idx
                        node_idx += 1
                    end
                end
            elseif startswith(line, "CTRIA3")
                if occursin(",", line)
                    parts = split(line, ",")
                    parts = strip.(parts)
                    if length(parts) >= 6
                        pid = parse(Int, parts[3])
                        g1 = parse(Int, parts[4])
                        g2 = parse(Int, parts[5])
                        g3 = parse(Int, parts[6])
                        
                        if haskey(grid_id_map, g1) && haskey(grid_id_map, g2) && haskey(grid_id_map, g3)
                            triangles[:, tri_idx] = [grid_id_map[g1], grid_id_map[g2], grid_id_map[g3]]
                            tri_tags[tri_idx] = pid
                            tri_idx += 1
                        end
                    end
                else
                    # Fixed width
                    pid_str = parse_nastran_field(line, 17, 24)
                    g1_str = parse_nastran_field(line, 25, 32)
                    g2_str = parse_nastran_field(line, 33, 40)
                    g3_str = parse_nastran_field(line, 41, 48)
                    
                    if !isempty(pid_str) && !isempty(g1_str) && !isempty(g2_str) && !isempty(g3_str)
                        pid = parse(Int, pid_str)
                        g1 = parse(Int, g1_str)
                        g2 = parse(Int, g2_str)
                        g3 = parse(Int, g3_str)
                        
                        if haskey(grid_id_map, g1) && haskey(grid_id_map, g2) && haskey(grid_id_map, g3)
                            triangles[:, tri_idx] = [grid_id_map[g1], grid_id_map[g2], grid_id_map[g3]]
                            tri_tags[tri_idx] = pid
                            tri_idx += 1
                        end
                    end
                end
            elseif startswith(line, "CTETRA")
                if occursin(",", line)
                    parts = split(line, ",")
                    parts = strip.(parts)
                    if length(parts) >= 7
                        # CTETRA, EID, PID, G1, G2, G3, G4
                        pid = parse(Int, parts[3])
                        g1 = parse(Int, parts[4])
                        g2 = parse(Int, parts[5])
                        g3 = parse(Int, parts[6])
                        g4 = parse(Int, parts[7])
                        
                        if haskey(grid_id_map, g1) && haskey(grid_id_map, g2) && haskey(grid_id_map, g3) && haskey(grid_id_map, g4)
                            tetras[:, tet_idx] = [grid_id_map[g1], grid_id_map[g2], grid_id_map[g3], grid_id_map[g4]]
                            tet_tags[tet_idx] = pid
                            tet_idx += 1
                        end
                    end
                else
                    # Fixed width
                    pid_str = parse_nastran_field(line, 17, 24)
                    g1_str = parse_nastran_field(line, 25, 32)
                    g2_str = parse_nastran_field(line, 33, 40)
                    g3_str = parse_nastran_field(line, 41, 48)
                    g4_str = parse_nastran_field(line, 49, 56)
                    
                    if !isempty(pid_str) && !isempty(g1_str) && !isempty(g2_str) && !isempty(g3_str) && !isempty(g4_str)
                        pid = parse(Int, pid_str)
                        g1 = parse(Int, g1_str)
                        g2 = parse(Int, g2_str)
                        g3 = parse(Int, g3_str)
                        g4 = parse(Int, g4_str)
                        
                        if haskey(grid_id_map, g1) && haskey(grid_id_map, g2) && haskey(grid_id_map, g3) && haskey(grid_id_map, g4)
                            tetras[:, tet_idx] = [grid_id_map[g1], grid_id_map[g2], grid_id_map[g3], grid_id_map[g4]]
                            tet_tags[tet_idx] = pid
                            tet_idx += 1
                        end
                    end
                end
            elseif startswith(line, "CHEXA")
                hex_nodes = _parse_chexa(line, f, grid_id_map)
                if !isnothing(hex_nodes)
                    hexas[:, hex_idx] = hex_nodes[1]
                    hex_tags[hex_idx] = hex_nodes[2]
                    hex_idx += 1
                end
            end
        end
        tri_idx, tet_idx, hex_idx
    end
    
    # Resize arrays if some elements were skipped
    if tri_idx <= trinum
        triangles = triangles[:, 1:tri_idx-1]
        tri_tags = tri_tags[1:tri_idx-1]
        trinum = tri_idx - 1
    end
    
    if tet_idx <= tetranum
        tetras = tetras[:, 1:tet_idx-1]
        tet_tags = tet_tags[1:tet_idx-1]
        tetranum = tet_idx - 1
    end

    if hex_idx <= hexanum
        hexas = hexas[:, 1:hex_idx-1]
        hex_tags = hex_tags[1:hex_idx-1]
        hexanum = hex_idx - 1
    end

    surface_mesh = TriangleMesh(trinum, node, triangles, tri_tags)
    volume_mesh = TetrahedraMesh(tetranum, node, tetras, tet_tags)
    hexa_mesh = HexahedraMesh(hexanum, node, hexas, hex_tags)
    
    return surface_mesh, volume_mesh, hexa_mesh
end

"""
    write_nas_mesh(pathname::String, mesh::TriangleMesh; scale=1.0)

Write a TriangleMesh to a Nastran (.nas) file.
"""
function write_nas_mesh(pathname::String, mesh::TriangleMesh; scale=1.0)
    open(pathname, "w") do f
        println(f, "\$ Created by EMSuite")
        println(f, "BEGIN BULK")
        
        # Write GRID points
        # GRID ID CP X Y Z CD PS SEID
        # Fixed width: 8 chars per field
        for i in 1:num_vertices(mesh)
            x = mesh.node[1, i] * scale
            y = mesh.node[2, i] * scale
            z = mesh.node[3, i] * scale
            
            # Using short format (free field with commas is easier and supported by read_nas_mesh)
            # But standard Nastran is fixed width. Let's use free field (comma separated) which is widely supported.
            # GRID, ID, , X, Y, Z
            println(f, "GRID, $i, , $x, $y, $z")
        end
        
        # Write CTRIA3 elements
        # CTRIA3 EID PID G1 G2 G3
        for i in 1:num_elements(mesh)
            n1 = mesh.triangles[1, i]
            n2 = mesh.triangles[2, i]
            n3 = mesh.triangles[3, i]
            pid = isempty(mesh.tags) ? 1 : mesh.tags[i]
            
            println(f, "CTRIA3, $i, $pid, $n1, $n2, $n3")
        end
        
        println(f, "ENDDATA")
    end
end
