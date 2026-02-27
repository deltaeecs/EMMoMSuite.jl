using ..CoreModule

function parse_nastran_field(line, start_col, end_col)
    if length(line) < start_col
        return ""
    end
    e = min(length(line), end_col)
    return strip(line[start_col:e])
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
    tri_tags = zeros(Int, trinum)
    tet_tags = zeros(Int, tetranum)
    
    # Second pass: read data
    tri_idx, tet_idx = open(pathname, "r") do f
        node_idx = 1
        tri_idx = 1
        tet_idx = 1
        
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
                        x = parse(FT, parts[4])
                        y = parse(FT, parts[5])
                        z = parse(FT, parts[6])
                        
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
                                x = parse(FT, x_str)
                                y = parse(FT, y_str)
                                z = parse(FT, z_str)
                                
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
                        x = parse(FT, x_str)
                        y = parse(FT, y_str)
                        z = parse(FT, z_str)
                        
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
            end
        end
        tri_idx, tet_idx
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

    if tetranum > 0
        return TetrahedraMesh(tetranum, node, tetras, tet_tags)
    else
        return TriangleMesh(trinum, node, triangles, tri_tags)
    end
end

"""
    read_mixed_nas_mesh(pathname::String; FT=Float64, scale=1.0)

Read a Nastran (.nas) mesh file and return both surface and volume meshes.
Returns `(surface_mesh::TriangleMesh, volume_mesh::TetrahedraMesh)`.
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
    tri_tags = zeros(Int, trinum)
    tet_tags = zeros(Int, tetranum)
    
    # Second pass: read data
    tri_idx, tet_idx = open(pathname, "r") do f
        node_idx = 1
        tri_idx = 1
        tet_idx = 1
        
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
                        x = parse(FT, parts[4])
                        y = parse(FT, parts[5])
                        z = parse(FT, parts[6])
                        
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
                                x = parse(FT, x_str)
                                y = parse(FT, y_str)
                                z = parse(FT, z_str)
                                
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
                        x = parse(FT, x_str)
                        y = parse(FT, y_str)
                        z = parse(FT, z_str)
                        
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
            end
        end
        tri_idx, tet_idx
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

    surface_mesh = TriangleMesh(trinum, node, triangles, tri_tags)
    volume_mesh = TetrahedraMesh(tetranum, node, tetras, tet_tags)
    
    return surface_mesh, volume_mesh
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
