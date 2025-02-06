```julia
using MakieStippling
using CUDA
#using Interpolations
using Random

using CairoMakie
using WGLMakie
#using MeshGrid
```

create an arbitrary density function

```julia
include("../simulate_densities.jl")
grid_densities = simulate_densities(timepoints = 10,grid_sz = (500,500),n_conditions=2)
```

```julia
trig = Observable(0)
 stipplemap(grid_densities,axis=(;aspect=1,limits=((0,500),(0,500)));one_more_iteration=trig)
# [trig[] = k for k in 1:10] # funny way to run 10 iterations
```

```julia


fig = Figure()
ax_dens = fig[1,1] = Axis(fig,title="underlying densities")
hidespines!(ax_dens)
hidedecorations!(ax_dens)
for k = 1:n_conditions
heatmap(fig[1,1][1,k],grid_densities[k])
hidedecorations!(current_axis())
end
fig
```

```julia
grid_single,sites_sets = MakieStippling.run_iterations!(grid_combined,sites,sites_grouping,grid_densities;n_iter=3,threshold_low=10e-10,threshold_high=0.65*10-6)


site_set_to_point2f(s::Set) = [Point2f(k[1],k[2]) for k in MakieStippling.linearized_sites.(collect(s),Ref(size(grid_combined)))]

fig = Figure()
ax = fig[1,1] = Axis(fig,aspect=1)
scatter!.(ax,site_set_to_point2f.(sites_sets);markerspace = :data,)
xlims!(ax,[0,500])
ylims!(ax,[0,500])
fig

```
