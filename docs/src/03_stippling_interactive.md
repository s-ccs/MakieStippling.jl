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



grid_combined, grid_singles, sites =
        setup(grid_densities, 100)

```

```julia
include("../simulate_densities.jl")
n_conditions = 3
grid_densities = CuArray.(simulate_densities(;timepoints = 10,grid_sz = (500,500),n_conditions))
colors = [:red,:green,:blue]
#---
#using TestImages
grid_densities = [CuArray(Float32.(1 .-rotr90(testimage("cameraman")))),CuArray(Float32.(rotr90(testimage("cameraman"))))]
colors = [:black,:white]

#---
_X = testimage("toucan")
using ImageTransformations
percentage_scale = 5
new_size = trunc.(Int, size(_X) .* percentage_scale)
X = rotr90(imresize(_X, new_size))

# exract a colorscheme
n = tempname()*".jpg";save(n,X)
colors = extract(n,5,25,0.1)
using Colors
function extract_color(X,color)
out = Colors.colordiff.(color,RGB.(X))

#out = out ./ maximum(out)
return Float32.(out)
end
colordiffs = extract_color.(Ref(X),colors)
m = maximum(maximum.(colordiffs))
X_cd = map(x->1 .-x./m,colordiffs)
X_cd = map(x-> (alpha.(X)).*x,X_cd)
grid_densities = CuArray.(X_cd)


#grid_densities = [CuArray(Float32.(rotr90(getproperty.(X,:r)))),
#                  CuArray(Float32.(rotr90(getproperty.(X,:g)))),
#                  CuArray(Float32.(rotr90(getproperty.(X,:b))))]
```

```julia
trig = Observable(20)
stipplemap(grid_densities;
    colors=colors[end:-1:1],
    density_factor=0.2,
    hysteresis=0.95,
    markersize=5,

    axis=(;aspect=1,limits=((0,size(grid_densities[1],1)),(0,size(grid_densities[1],2)))),
    one_more_iteration=trig)
```

run setup without iterations

```julia

grid_combined, grid_singles, sites =
        MakieStippling.setup(grid_densities, 100)
MakieStippling.jfa_voronoi_gpu!(grid_combined, collect(reduce(union, to_value.(sites))))
n_groups = length(grid_densities)
[MakieStippling.jfa_voronoi_gpu!(grid_singles[k], collect(to_value(sites[k]))) for k = 1:n_groups]


_sites = Set.(MakieStippling.linearized_sites.(collect.(to_value.(sites)), Ref(size(grid_combined))))

@run MakieStippling.one_iteration!(
            grid_combined,
            grid_singles,
            _sites,
            to_value(grid_densities);
            markersize = 5,
            hysteresis_a = 0.8
        )


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
