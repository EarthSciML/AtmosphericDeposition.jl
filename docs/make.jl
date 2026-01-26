using AtmosphericDeposition
using Documenter

DocMeta.setdocmeta!(
    AtmosphericDeposition,
    :DocTestSetup,
    :(using AtmosphericDeposition);
    recursive = true
)

makedocs(;
    modules = [AtmosphericDeposition],
    authors = "EarthSciML authors and contributors",
    repo = "https://github.com/EarthSciML/AtmosphericDeposition.jl/blob/{commit}{path}#{line}",
    sitename = "AtmosphericDeposition.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://deposition.earthsci.dev",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Dry Deposition" => [
            "Overview" => "dry_deposition.md",
            "Wesley" => "wesley1989.md",
            "Water Surfaces" => "water_deposition.md"],
        "Wet Deposition" => ["EMEP" => "emep.md"],
        "API" => "api.md"
    ]
)

deploydocs(; repo = "github.com/EarthSciML/AtmosphericDeposition.jl")
