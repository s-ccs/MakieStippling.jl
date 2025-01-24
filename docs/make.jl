using MakieStippling
using Documenter

DocMeta.setdocmeta!(MakieStippling, :DocTestSetup, :(using MakieStippling); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [MakieStippling],
    authors = "Benedikt Ehinger <benedikt.ehinger@vis.uni-stuttgart.de>",
    repo = "https://github.com/s-ccs/MakieStippling.jl/blob/{commit}{path}#{line}",
    sitename = "MakieStippling.jl",
    format = Documenter.HTML(; canonical = "https://s-ccs.github.io/MakieStippling.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/s-ccs/MakieStippling.jl")
