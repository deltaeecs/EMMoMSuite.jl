using Markdown

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

files = String[]
for (root, _, names) in walkdir(DOCS_SRC)
    for name in names
        endswith(name, ".md") || continue
        push!(files, joinpath(root, name))
    end
end

for file in sort(files)
    text = read(file, String)
    converted = convert_markdown_for_documenter(text)
    try
        Markdown.parse(converted)
        println("OK  " * relpath(file, DOCS_SRC))
    catch err
        println("FAIL " * relpath(file, DOCS_SRC))
        showerror(stdout, err)
        println()
        rethrow(err)
    end
end

println("ALL_PARSED")