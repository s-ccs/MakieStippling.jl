module MakieStippling
using Distributions: value_support

using Interpolations, Random
using Distributions

using Distances
using KernelAbstractions

using Atomix

using Makie

include("jfa_voronoi.jl")
include("grid_features.jl")
include("density.jl")
include("stippling.jl")
include("stipplemap.jl")

export stipplemap, stipplemap!
end
