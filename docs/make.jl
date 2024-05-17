using NICE
using Documenter

DocMeta.setdocmeta!(NICE, :DocTestSetup, :(using NICE); recursive=true)

makedocs(;
    modules=[NICE],
    authors="QC-Devs",
    repo="https://github.com/theochem/NICE.jl/blob/{commit}{path}#{line}",
    sitename="NICE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nice.qcdevs.org/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/theochem/NICE.jl",
    devbranch="main",
)
