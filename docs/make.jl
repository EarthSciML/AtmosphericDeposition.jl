using DepositionMTK
using Documenter

DocMeta.setdocmeta!(DepositionMTK, :DocTestSetup, :(using DepositionMTK); recursive=true)

makedocs(;
    modules=[DepositionMTK],
    authors="Chris Tessum <ctessum@gmail.com> and contributors",
    repo="https://github.com/ctessum/DepositionMTK.jl/blob/{commit}{path}#{line}",
    sitename="DepositionMTK.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ctessum.github.io/DepositionMTK.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ctessum/DepositionMTK.jl",
)
