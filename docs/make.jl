using ElectrochemicalReactor
using Documenter

DocMeta.setdocmeta!(ElectrochemicalReactor, :DocTestSetup, :(using ElectrochemicalReactor); recursive=true)

makedocs(;
    modules=[ElectrochemicalReactor],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/ElectrochemicalReactor.jl/blob/{commit}{path}#{line}",
    sitename="ElectrochemicalReactor.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/ElectrochemicalReactor.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/ElectrochemicalReactor.jl",
    devbranch="main",
)
