using AtmosphericDeposition
using AtmosphericDeposition: DryDeposition
using Documenter

DocMeta.setdocmeta!(
    AtmosphericDeposition,
    :DocTestSetup,
    :(using AtmosphericDeposition; using AtmosphericDeposition.DryDeposition);
    recursive = true
)

makedocs(;
    modules = [AtmosphericDeposition, DryDeposition],
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
            "Seinfeld & Pandis Ch. 19" => "seinfeld_pandis_ch19.md",
        ],
        "Wet Deposition" => [
            "EMEP" => "emep.md",
            "Seinfeld & Pandis (2006) Ch. 20" => "seinfeld_pandis_2006_ch20.md",
            "Luo & Yu (2023)" => "luo2023.md",
        ],
        "API" => "api.md",
    ]
)

deploydocs(; repo = "github.com/EarthSciML/AtmosphericDeposition.jl")
