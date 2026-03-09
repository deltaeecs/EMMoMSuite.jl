using Documenter
using EMSuite

makedocs(
    sitename = "EMSuite.jl",
    format = Documenter.HTML(repolink = ""),
    modules = [EMSuite],
    remotes = nothing,
    checkdocs = :none,      # 不强制要求所有 docstring 出现在手册中
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
    ]
)

# deploydocs(
#     repo = "github.com/username/EMSuite.jl.git",
# )
