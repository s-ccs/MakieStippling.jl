```julia
using MakieStippling
using CUDA
using DiscreteVoronoi
using StaticArrays
using CairoMakie


```

```julia
grid = Int.(zeros(30,30))
grid_static = zeros(SVector{2,Int}, size(grid,1),size(grid,2))
cu_grid = CuArray(grid)

n_sites = 100
sites = CartesianIndices(size(grid))[rand(1:length(grid),n_sites)]
#sites_lin = MakieStippling.linearized_sites(sites, size(grid))

sites_static = [SVector{2,Int}(x,y) for (x,y) in Tuple.(sites)]
#sites = [CartesianIndex(1,1),CartesianIndex(10,10)]
#grid[sites] .= 1:size(sites,1)
#---
MakieStippling.jfa_voronoi_gpu!(grid,sites)
MakieStippling.jfa_voronoi_gpu!(cu_grid,sites)

DiscreteVoronoi.jfa_voronoi!(grid_static,sites_static)

# convert to linearized
grid_static_lin =MakieStippling.linearized_sites(grid_static,size(grid_static))

n_elems = *(size(grid)...)

# DiscreteVoronoi vs. CPU
 sum(grid_static_lin .== grid) / n_elems
```

```julia
# CPU vs GPU
sum(Matrix(cu_grid) .== grid) / n_elems
```

```julia
# DiscreteVoronoi vs. GPU
sum(Matrix(cu_grid) .== grid_static_lin)/n_elems

```

```julia
f = Figure()

heatmap(f[1,1],grid_static_lin,colormap=:tab10,axis=(;title="DiscreteVoronoi.jl"))
heatmap(f[1,2],grid,colormap=:tab10,axis=(;title="cpu"))
heatmap(f[1,3],Matrix(cu_grid),colormap=:tab10,axis=(;title="cu"))

heatmap(f[2,1],grid.==grid_static_lin,axis=(;title="grid - DiscreteVoronoi.jl"))
heatmap(f[2,2],Matrix(cu_grid).==grid_static_lin,axis=(;title="Cu - DiscreteVoronoi.jl"))
heatmap(f[2,3],Matrix(cu_grid).==grid,axis=(;title="Cu - Grid"))

f
