using MaximumEntropyMomentClosures
using Documenter

DocMeta.setdocmeta!(MaximumEntropyMomentClosures, :DocTestSetup, :(using MaximumEntropyMomentClosures); recursive=true)

makedocs(;
    modules=[MaximumEntropyMomentClosures],
    authors="roman-schaerer <r.p.s@posteo.ch> and contributors",
    repo="https://github.com/roman-schaerer/MaximumEntropyMomentClosures.jl/blob/{commit}{path}#{line}",
    sitename="MaximumEntropyMomentClosures.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://roman-schaerer.github.io/MaximumEntropyMomentClosures.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/roman-schaerer/MaximumEntropyMomentClosures.jl",
    devbranch="main",
)
