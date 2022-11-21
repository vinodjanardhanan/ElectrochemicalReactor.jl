using ElectrochemicalCell
using Documenter

DocMeta.setdocmeta!(ElectrochemicalCell, :DocTestSetup, :(using ElectrochemicalCell); recursive=true)

makedocs(;
    modules=[ElectrochemicalCell],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/ElectrochemicalCell.jl/blob/{commit}{path}#{line}",
    sitename="ElectrochemicalCell.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/ElectrochemicalCell.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/ElectrochemicalCell.jl",
    devbranch="main",
)
