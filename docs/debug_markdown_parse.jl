using Markdown

include("make.jl")

src_dir = joinpath(@__DIR__, "src")
converted_dir = prepare_documenter_source(src_dir)

failures = String[]
for (root, _, files) in walkdir(converted_dir)
    for file in files
        endswith(file, ".md") || continue
        path = joinpath(root, file)
        content = read(path, String)
        try
            Markdown.parse(content)
        catch err
            push!(failures, path)
            println("FAILED: ", path)
            showerror(stdout, err)
            println()
        end
    end
end

println("FAILURE_COUNT=", length(failures))
for path in failures
    println(path)
end
