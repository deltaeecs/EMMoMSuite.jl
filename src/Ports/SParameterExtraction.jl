"""
    SParameterExtraction — S 参数提取与 Touchstone I/O

# S 矩阵提取（从 MoM 阻抗矩阵）

采用端口 Z 矩阵法：从全阻抗矩阵 Z_efie 提取端口索引处的子矩阵 Z_p，
然后通过标准公式转换为 S 矩阵：

    S = (Z_p − Z_ref · I) · (Z_p + Z_ref · I)⁻¹

对多端口情况，Z_p 是 N_port×N_port 子矩阵，支持不同端口参考阻抗。

# Touchstone 格式（v1.0）

```
!注释
# Hz S RI R 50
freq  S11_real S11_imag  S21_real S21_imag  ...
```

行顺序：端口 (j,k) 以 j 快、k 慢的方式排列（Touchstone 标准），
即对于 2 端口：S11 S21 S12 S22（每列先），S 分量顺序为列主序。

# 参考
Pozar《Microwave Engineering》附录 A；Touchstone® File Format Specification v1.1。
"""

# ─────────────────────────────────────────────────────────────────────────────
# SParameterData 数据结构
# ─────────────────────────────────────────────────────────────────────────────

"""
    SParameterData

多端口 S 参数数据容器。

# 字段
- `frequencies`:      频率向量（Hz）
- `s_matrices`:       每个频率点的 S 矩阵（N×N ComplexF64）列表
- `port_impedances`:  每个端口的参考阻抗（Ω；标量表示所有端口相同）
- `comments`:         注释字符串列表（对应 Touchstone '!' 注释行）
"""
struct SParameterData
    frequencies::Vector{Float64}
    s_matrices::Vector{Matrix{ComplexF64}}
    port_impedances::Vector{Float64}
    comments::Vector{String}
end

"""
    sparam_data(freqs, smats, port_z, cmts) → SParameterData

带验证的 SParameterData 工厂函数。自动将输入转换为正确的元素类型。
"""
function sparam_data(freqs, smats, port_z, cmts = String[])
    n_f = length(freqs)
    sv  = collect(Matrix{ComplexF64}, smats)
    @assert length(sv) == n_f "s_matrices 数量须与 frequencies 匹配"
    if !isempty(sv)
        n_p = size(sv[1], 1)
        @assert size(sv[1], 2) == n_p "S 矩阵必须是方阵"
        @assert length(port_z) == n_p "port_impedances 长度须等于端口数"
    end
    return SParameterData(
        collect(Float64, freqs),
        sv,
        collect(Float64, port_z),
        collect(String, cmts),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# S 矩阵提取
# ─────────────────────────────────────────────────────────────────────────────

"""
    extract_s_matrix(ports, Z_efie; Z_ref=50.0) → Matrix{ComplexF64}

从 MoM 全阻抗矩阵提取端口 S 矩阵。

# 算法
1. 提取端口子矩阵：Z_p[j,k] = Z_efie[ports[j].edge_idx, ports[k].edge_idx]
2. 转换为 S 矩阵：S = (Z_p − Z_ref·I)·(Z_p + Z_ref·I)⁻¹

# 精度说明
当端口 DOF 之间存在非端口 DOF 耦合时，应使用 Schur 补消去法求得精确 Z_p。
当前实现为直接子矩阵提取（适用于端口间距 >> λ 的情况）。

# 参数
- `ports`:   `LumpedPort` 列表（指定端口边索引）
- `Z_efie`:  全阻抗矩阵（N_basis × N_basis）
- `Z_ref`:   参考阻抗（Ω，默认 50.0；所有端口相同）
"""
function extract_s_matrix(
    ports::AbstractVector{<:AbstractPort},
    Z_efie::AbstractMatrix{<:Complex};
    Z_ref::Real = 50.0,
)
    # 提取每个 LumpedPort 的边索引
    edge_idx = [_port_edge_idx(p) for p in ports]
    N = length(edge_idx)

    # 端口 Z 矩阵（子矩阵提取）
    Z_p = Z_efie[edge_idx, edge_idx]

    # Z → S 转换：S = (Z − Z0) · (Z + Z0)⁻¹
    # 使用 A/B（即 A·B⁻¹）避免显式 inv() 以提高数值稳定性
    Z0I = Z_ref * I
    S = (Z_p - Z0I) / (Z_p + Z0I)
    return S
end

"""
    _port_edge_idx(port) → Int

从端口对象提取边索引。支持 LumpedPort；对其他端口类型抛出友好错误。
"""
function _port_edge_idx(p::LumpedPort)
    return p.edge_idx
end
_port_edge_idx(p::AbstractPort) =
    error("extract_s_matrix: port type $(typeof(p)) 尚不支持直接 S 矩阵提取")

# ─────────────────────────────────────────────────────────────────────────────
# Touchstone 1.0 写入
# ─────────────────────────────────────────────────────────────────────────────

"""
    write_touchstone(data::SParameterData, path::String; version::Int=1)

将 SParameterData 写出为 Touchstone 1.0 格式文件（.s1p/.s2p/.sNp）。

# 格式
```
!注释（来自 data.comments）
# Hz S RI R <Z_ref>
<freq_Hz>  <S11r S11i  S21r S21i  S12r S12i  S22r S22i>  （列主序）
```

# 端口顺序（Touchstone 标准）
对于 N 端口，每行写出所有 S 分量，顺序为列主序：
S11 S21 ... SN1  S12 S22 ... SN2  ... SN1 ... SNN
（即先列 j=1 的所有 i，再列 j=2 的所有 i...）

但 2 端口 Touchstone 1.0 的惯例：每行写 S11 S21 S12 S22，
即 S[j,k] 按 (1,1) (2,1) (1,2) (2,2) 顺序交叉排列。

这里实现 **标准 Touchstone 1.0** 行格式：
- 1 端口：S11
- 2 端口：S11 S21 S12 S22
- 通用 N 端口：S11...SN1 S12...SN2 ... S1N...SNN（列主序）

注：部分 EDA 工具对 2 端口的行/列主序解释不同，
上述实现与 Keysight 的 Touchstone 1.1 规范保持一致。
"""
function write_touchstone(data::SParameterData, path::String; version::Int = 1)
    isempty(data.port_impedances) && error("write_touchstone: port_impedances 不能为空")
    Z_ref  = data.port_impedances[1]   # 假设所有端口相同 Z_ref
    n_port = isempty(data.s_matrices) ? length(data.port_impedances) : size(data.s_matrices[1], 1)

    open(path, "w") do io
        # 注释行
        for c in data.comments
            println(io, c)
        end
        # 选项行
        println(io, "# Hz S RI R $(Z_ref)")
        # 数据行
        for (freq, S) in zip(data.frequencies, data.s_matrices)
            # 列主序：先列出 S[:,1]，再 S[:,2]，...
            tokens = Float64[]
            for j in 1:n_port
                for i in 1:n_port
                    push!(tokens, real(S[i, j]))
                    push!(tokens, imag(S[i, j]))
                end
            end
            print(io, freq)
            for t in tokens
                print(io, "  ", t)
            end
            println(io)
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Touchstone 1.0 读取
# ─────────────────────────────────────────────────────────────────────────────

"""
    read_touchstone(path::String) → SParameterData

读取 Touchstone 1.0 格式文件，返回 SParameterData。

支持：
- `!` 注释行
- `# Hz S RI R <Z_ref>` 选项行（当前仅支持 Hz、S 参数、RI 格式）
- 多端口数据行（列主序，与 write_touchstone 对应）

限制：
- 仅支持 RI（实虚）格式；MA（幅相）格式暂不支持
- 频率单位必须为 Hz
"""
function read_touchstone(path::String)
    comments     = String[]
    freq_list    = Float64[]
    smats        = Matrix{ComplexF64}[]
    Z_ref        = 50.0
    n_port       = 0   # 从文件扩展名推断

    # 从文件扩展名推断端口数
    ext = lowercase(splitext(path)[2])   # ".s2p" → ".s2p"
    if startswith(ext, ".s") && endswith(ext, "p")
        n_str = ext[3:end-1]
        n_port = tryparse(Int, n_str)
        isnothing(n_port) && (n_port = 0)
    end

    open(path, "r") do io
        for line in eachline(io)
            line_stripped = strip(line)
            isempty(line_stripped) && continue

            # 注释行
            if startswith(line_stripped, '!')
                push!(comments, line_stripped)
                continue
            end

            # 选项行（以 '#' 开头）
            if startswith(line_stripped, '#')
                tokens = split(lowercase(line_stripped))
                # # Hz S RI R <Z_ref>
                for (ti, tok) in enumerate(tokens)
                    if tok == "r" && ti < length(tokens)
                        Z_ref = parse(Float64, tokens[ti+1])
                    end
                end
                # 仅支持 Hz + RI
                # （MA、DB 格式验证暂时跳过，不报错）
                continue
            end

            # 数据行：第一个数为频率，后续为 S 参数（实部 虚部 交替）
            tokens = split(line_stripped)
            isempty(tokens) && continue

            freq = parse(Float64, tokens[1])
            push!(freq_list, freq)

            # 解析 S 参数对
            vals  = [parse(Float64, t) for t in tokens[2:end]]
            n_val = length(vals)
            n_pairs = n_val ÷ 2
            if n_port == 0
                n_port = isqrt(n_pairs)  # 从数据推断 n×n 矩阵
                n_port * n_port != n_pairs &&
                    @warn "Touchstone 行数据含 $n_pairs 对值，非完全平方数，端口数推断为 $n_port（可能不正确）"
            end

            N  = n_port
            S  = zeros(ComplexF64, N, N)
            k  = 1
            for j in 1:N
                for i in 1:N
                    if k + 1 <= n_val
                        S[i, j] = complex(vals[k], vals[k+1])
                        k += 2
                    end
                end
            end
            push!(smats, S)
        end
    end

    port_z = fill(Z_ref, n_port)
    return SParameterData(freq_list, smats, port_z, comments)
end
