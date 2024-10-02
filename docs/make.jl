using VectorFitting
using Documenter

DocMeta.setdocmeta!(VectorFitting, :DocTestSetup, :(using VectorFitting); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [VectorFitting],
    authors = "Pedro H. N. Vieira <phnvieira@proton.me>",
    repo = "https://github.com/pedrohnv/VectorFitting.jl/blob/{commit}{path}#{line}",
    sitename = "VectorFitting.jl",
    format = Documenter.HTML(; canonical = "https://pedrohnv.github.io/VectorFitting.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/pedrohnv/VectorFitting.jl")
