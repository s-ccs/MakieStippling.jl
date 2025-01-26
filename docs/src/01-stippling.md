```julia
using MakieStippling
using CUDA
#using Interpolations
using Random

using CairoMakie
using MeshGrid
```

create random data

```julia
n_sites = 5000
n_conditions = 2
sites_grouping = rand(1:n_conditions,n_sites)


grid = Int.(zeros(500,500))
grid_combined = CuArray(grid)
sites = CartesianIndices(size(grid))[rand(1:length(grid),n_sites)]
nothing #hide
```

create an arbitrary density function

```julia
times = collect(range(0,size(grid_combined,1),length=5))
density_interpolators = Vector{Any}(undef,n_conditions)
max_density = 1
for k = 1:n_conditions
    mean_data = rand(MersenneTwister(k),5).*size(grid,1)/4 .+ size(grid,1)/2
    std_data = [20,10,5,15,30.]#abs.(rand(MersenneTwister(k*2),5)).*size(grid,1)

    itp_mean = MakieStippling.Interpolations.interpolate((times,), mean_data, MakieStippling.Interpolations.Gridded(MakieStippling.Interpolations.Linear()))
    itp_std = MakieStippling.Interpolations.interpolate((times,), std_data, MakieStippling.Interpolations.Gridded(MakieStippling.Interpolations.Linear()))

    density_interpolators[k] = MakieStippling.DensityInterpolator(itp_mean,itp_std,max_density)

end

grid_densities = Vector(undef,n_conditions)
for k = 1:n_conditions
grid_densities[k] = Array{Float64}(undef,size(grid_combined))
for x = 1:size(grid_density,1)
for y = 1:size(grid_density,2)

grid_densities[k][x,y] =  MakieStippling.calculate_density(density_interpolators[k], (x,y))
end
end
end


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
