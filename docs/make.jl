using Documenter
using AstroCoords

makedocs(;
    modules=[AstroCoords],
    format=Documenter.HTML(;
        prettyurls=!("local" in ARGS), highlights=["yaml"], ansicolor=true
    ),
    sitename="AstroCoordss.jl",
    authors="Jordan Murphy",
    pages=[
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "API" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(; repo="github.com/jmurphy6895/AstroCoords.jl.git", target="build")
