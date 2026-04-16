using Documenter
using TriShellFiniteElement

makedocs(
    sitename = "TriShellFiniteElement.jl",
    modules = [TriShellFiniteElement],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/runtosolve/TriShellFiniteElement.jl.git",
)
