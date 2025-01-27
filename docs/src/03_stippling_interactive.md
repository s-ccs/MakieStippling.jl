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
grid_densities = simulate_densities(timepoints = 5,grid_sz = size(grid_combined),n_conditions=2)

n_conditions = 1
grid_densities = [testimage("cameraman")]
```

run setup without iterations

```julia
grid_single,sites_sets = MakieStippling.run_iterations!(grid_combined,sites,sites_grouping,grid_densities;
        n_iter=0,threshold_low=10e-7,threshold_high=10e-6)
```

setup `Observables.jl` to quickly update the plot

```julia
site_set_to_point2f(s::Set) = [Point2f(k[1],k[2]) for k in MakieStippling.linearized_sites.(collect(s),Ref(size(grid_combined)))]

sites_sets_obs = Observable(sites_sets)
sites_sets_obs_point = [@lift(site_set_to_point2f($(sites_sets_obs)[k])) for k = 1:n_conditions]

sites_combined_obs = @lift(reduce(union,$(sites_sets_obs)))
grid_combined_obs = Observable(Matrix(grid_combined))
```

add a function that does one step with a given threshold

```julia
function update_observables!(grid_combined_obs,sites_sets_obs,grid_combined,grid_single,sites_sets,grid_densities; threshold_low=10e-15,threshold_high=0.0000005)
MakieStippling.one_iteration!(grid_combined,grid_single,sites_sets,grid_densities;
        threshold_low,threshold_high)

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
scatter!.(ax_n,sites_sets_obs_point,markerspace = :data)
sg = SliderGrid(
    fig_slider[2,1][1,1],
    (label = L"thres_{low}", range =0.001:0.001:1, format = "{:.3f}", startvalue = 0.01),
    (label = L"thres_{high}", range =0.001:0.001:1, format = "{:.3f}", startvalue = 0.1),
    tellwidth=false,)
    #height=100)
    #tellheight = false,tellwidth=false)


denseXareas_obs = map(sites_sets_obs) do s
    _tmp =MakieStippling.grid_features(grid_singles,reduce(union,s))
    _tmp[1] .* grid_densities[1]
    return _tmp[1] # areas
end

axh = hist(fig_slider[2,1][1,2],Float64.(grid_densities[1][:]))
#axh.axis.xticks = []


reset_button =fig_slider[3,1] = WGLMakie.Makie.Button(fig_slider, label = "restart",tellwidth=false)
on(reset_button.clicks) do n
    println("resetting")
    grid_single_new,sites_sets_new = MakieStippling.run_iterations!(grid_combined,sites,sites_grouping,grid_densities;
    n_iter=0,threshold_low=10e-7,threshold_high=-0.5*10e-6)
    grid_single .= grid_single_new
    sites_sets .= sites_sets_new
    set_close_to!(sg.sliders[1], -7)
    set_close_to!(sg.sliders[2], -6.5)
    update_observables!(grid_combined_obs,sites_sets_obs,grid_combined,grid_single,sites_sets,grid_densities)
end

iterate_button =fig_slider[3,1][1,2] = WGLMakie.Makie.Button(fig_slider, label = "one iteration",tellwidth=false)
on(iterate_button.clicks) do n
    update_observables!(grid_combined_obs,sites_sets_obs,grid_combined,grid_single,sites_sets,grid_densities,threshold_low=10^sg.sliders[1].value.val,threshold_high=10^sg.sliders[2].value.val)
end
fig_slider
```
