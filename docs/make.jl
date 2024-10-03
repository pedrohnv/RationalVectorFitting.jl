using RationalVectorFitting
using Documenter

DocMeta.setdocmeta!(
    RationalVectorFitting,
    :DocTestSetup,
    :(using RationalVectorFitting);
    recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [RationalVectorFitting],
    authors = "Pedro H. N. Vieira <phnvieira@proton.me>",
    repo = "https://github.com/pedrohnv/RationalVectorFitting.jl/blob/{commit}{path}#{line}",
    sitename = "RationalVectorFitting.jl",
    format = Documenter.HTML(;
        canonical = "https://pedrohnv.github.io/RationalVectorFitting.jl",
    ),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/pedrohnv/RationalVectorFitting.jl")
