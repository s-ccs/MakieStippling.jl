# Interactive example

this potentially does not run via the CLI because `WGLMakie` might not work. I haven't checked.

```julia
using MakieStippling
using CUDA
#using Interpolations
using Random
using WGLMakie
using MeshGrid
```

```julia
n_sites = 5000
n_conditions = 2
sites_grouping = rand(1:n_conditions,n_sites)


grid = Int.(zeros(500,500))
grid_combined = CuArray(grid)
sites = CartesianIndices(size(grid))[rand(1:length(grid),n_sites)]
```

```julia
include("../simulate_densities.jl")
n_conditions = 2
grid_densities = simulate_densities(;timepoints = 10,grid_sz = size(grid_combined),n_conditions)

#grid_densities = [testimage("cameraman")]
```

run setup without iterations

```julia
grid_singles,sites_sets = MakieStippling.run_iterations!(grid_combined,sites,sites_grouping,grid_densities;
        n_iter=0)
```

setup `Observables.jl` to quickly update the plot

```julia
site_set_to_point2f(s::Set) = [Point2f(k[1],k[2]) for k in MakieStippling.linearized_sites.(collect(s),Ref(size(grid_combined)))]

sites_sets_obs = Observable(sites_sets)
sites_sets_obs_point = [@lift(site_set_to_point2f($(sites_sets_obs)[k])) for k = 1:n_conditions]

sites_combined_obs = @lift(reduce(union,$(sites_sets_obs)))
grid_combined_obs = Observable(Matrix(grid_combined))
markersize_obs = Observable(5)
```

add a function that does one step with a given threshold

```julia
function update_observables!(grid_combined_obs,sites_sets_obs,grid_combined,grid_singles,sites_sets,grid_densities; hysteresis_a,markersize)
MakieStippling.one_iteration!(grid_combined,grid_singles,sites_sets,grid_densities;markersize,
        hysteresis_a)

    grid_combined_obs.val = Matrix(grid_combined)
    sites_sets_obs.val= sites_sets
    notify(sites_sets_obs)
    notify(grid_combined_obs)
end
```

setup the figure

```julia
fig_slider = Figure(size=(800,400))
ax_n = fig_slider[1,1] = Axis(fig_slider,aspect=1)
scatter!.(ax_n,sites_sets_obs_point,markerspace = :data,markersize=markersize_obs)
sg = SliderGrid(
    fig_slider[2,1][1,1],
    (label = L"hysteresis_a", range =0:0.05:1, format = "{:.2f}", startvalue = 0.6),
    (label = L"markersize", range =1:1:20, format = "{:d}", startvalue = 10),

    tellwidth=false,)
    #height=100)
    #tellheight = false,tellwidth=false)

on(sg.sliders[2].value) do s
    markersize_obs[] = s
    end

#denseXareas_obs = map(sites_sets_obs) do s
#    _tmp =MakieStippling.grid_features(grid_singles[1],s[1])
#    _tmp[1] .* grid_densities[1]
#    return _tmp[1] # areas
#end

#axh = hist(fig_slider[2,1][1,2],Float64.(grid_densities[1][:]))
#axh.axis.xticks = []


reset_button =fig_slider[3,1] = WGLMakie.Makie.Button(fig_slider, label = "restart",tellwidth=false)
on(reset_button.clicks) do n
    println("resetting")
    grid_singles_new,sites_sets_new = MakieStippling.run_iterations!(grid_combined,sites,sites_grouping,
    grid_densities;
    n_iter=0,markersize=10,hysteresis_a=0.6)
    grid_singles .= grid_singles_new
    sites_sets .= sites_sets_new
    set_close_to!(sg.sliders[1], 0.6)
    set_close_to!(sg.sliders[2], 10)
    #set_close_to!(sg.sliders[2], -6.5)
    update_observables!(grid_combined_obs,
    sites_sets_obs,grid_combined,grid_singles,sites_sets,grid_densities;
    markersize = to_value(markersize_obs),
    hysteresis_a=0.6)
end

iterate_button =fig_slider[3,1][1,2] = WGLMakie.Makie.Button(fig_slider, label = "one iteration",tellwidth=false)
on(iterate_button.clicks) do n
    update_observables!(grid_combined_obs,sites_sets_obs,
    grid_combined,grid_singles,sites_sets,
    grid_densities;
    markersize = to_value(markersize_obs),
    hysteresis_a=sg.sliders[1].value.val)
end


[heatmap(fig_slider[1,2][1,k],Float64.(grid_densities[k])) for k in 1:length(grid_singles)]
fig_slider
```
