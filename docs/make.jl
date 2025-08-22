using Documenter
using FrequencyMaxwell

makedocs(;
    modules=[FrequencyMaxwell],
    authors="Dohyeon Lee <dleh428@kaist.ac.kr>",
    repo="https://github.com/your-username/FrequencyMaxwell.jl/blob/{commit}{path}#{line}",
    sitename="FrequencyMaxwell.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://your-username.github.io/FrequencyMaxwell.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "manual/getting_started.md",
            "Configuration" => "manual/configuration.md",
            "Solvers" => "manual/solvers.md",
            "Sources" => "manual/sources.md",
            "Phantoms" => "manual/phantoms.md",
        ],
        "Examples" => [
            "Basic Scattering" => "examples/basic_scattering.md",
            "Phantom Gallery" => "examples/phantom_gallery.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/your-username/FrequencyMaxwell.jl",
    devbranch="main",
)
