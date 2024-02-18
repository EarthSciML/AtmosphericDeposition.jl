using AtmosphericDeposition
using Documenter

DocMeta.setdocmeta!(AtmosphericDeposition, :DocTestSetup, :(using AtmosphericDeposition); recursive=true)

makedocs(;
    modules=[AtmosphericDeposition],
    authors="EarthSciML authors and contributors",
    repo="https://github.com/earthsciml/AtmosphericDeposition.jl/blob/{commit}{path}#{line}",
    sitename="AtmosphericDeposition.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://earthciml.github.io/AtmosphericDeposition.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/earthsciml/AtmosphericDeposition.jl",
)
