using NICE
using Documenter

DocMeta.setdocmeta!(NICE, :DocTestSetup, :(using NICE); recursive=true)

makedocs(;
    modules=[NICE],
    authors="QC-Devs",
    repo="https://github.com/quantumelephant/NICE.jl/blob/{commit}{path}#{line}",
    sitename="NICE.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://quantumelephant.github.io/NICE.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/quantumelephant/NICE.jl",
    devbranch="main",
)
