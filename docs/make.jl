using Documenter
using EMMoMSuite

const DOCS_ROOT = @__DIR__
const DOCS_SRC = joinpath(DOCS_ROOT, "src")

function convert_inline_math(line::AbstractString)
    replace(line, r"(?<!\\)\$(?!\$)(.+?)(?<!\\)\$(?!\$)" => s"``\1``")
end

function convert_markdown_for_documenter(text::AbstractString)
    lines = split(text, '\n'; keepempty = true)
    output = String[]
    math_lines = String[]
    in_fence = false
    in_display_math = false

    for line in lines
        stripped = strip(line)

        if in_display_math
            if stripped == "\$\$"
                push!(output, "```math")
                append!(output, math_lines)
                push!(output, "```")
                empty!(math_lines)
                in_display_math = false
            else
                push!(math_lines, line)
            end
            continue
        end

        if startswith(stripped, "```")
            push!(output, line)
            in_fence = !in_fence
            continue
        end

        if !in_fence && stripped == "\$\$"
            in_display_math = true
            continue
        end

        push!(output, in_fence ? line : convert_inline_math(line))
    end

    join(output, '\n')
end

function prepare_documenter_source(src_dir::AbstractString)
    tmp_root = mktempdir()
    out_dir = joinpath(tmp_root, "src")
    mkpath(out_dir)

    for (root, _, files) in walkdir(src_dir)
        rel = relpath(root, src_dir)
        dst_root = rel == "." ? out_dir : joinpath(out_dir, rel)
        mkpath(dst_root)

        for file in files
            src_path = joinpath(root, file)
            dst_path = joinpath(dst_root, file)
            if endswith(file, ".md")
                content = read(src_path, String)
                write(dst_path, convert_markdown_for_documenter(content))
            else
                cp(src_path, dst_path; force = true)
            end
        end
    end

    out_dir
end

generated_source = prepare_documenter_source(DOCS_SRC)

makedocs(
    sitename = "EMMoMSuite.jl",
    source = generated_source,
    format = Documenter.HTML(repolink = "", prettyurls = false, edit_link = "master"),
    modules = [EMMoMSuite],
    remotes = nothing,
    checkdocs = :none,      # Do not require every docstring to appear in the manual.
    warnonly = true,
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Installation" => "guide/installation.md",
            "Quick Start" => "guide/quick_start.md",
            "Advanced" => "guide/advanced.md",
            "Examples" => "guide/examples.md",
        ],
        "Theory" => [
            "Overview" => "theory/overview.md",
            "Electromagnetics" => "theory/electromagnetics.md",
            "Integral Equations" => "theory/integral_equations.md",
            "Basis Functions" => "theory/basis_functions.md",
            "Method of Moments" => "theory/method_of_moments.md",
            "Fast Algorithms" => "theory/fast_algorithms.md",
            "Excitations & Ports" => "theory/excitations.md",
            "Solvers" => "theory/solvers.md",
            "Post-Processing" => "theory/post_processing.md",
        ],
        "API Reference" => "api/public_api.md",
    ],
)

# deploydocs(
#     repo = "github.com/username/EMMoMSuite.jl.git",
# )
