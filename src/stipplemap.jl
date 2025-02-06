@recipe(StippleMap, densities) do scene
    Attributes(
        markersize = 10.0,
        hysteresis = 0.6,
        density_factor = 0.5,
        n_initial_sites = 1000,
        n_iter = 10,
        one_more_iteration = 0,
        colors = nothing,
    )
end

function setup(grid_densities, initial_sites)
    n_groups = length(grid_densities)
    grid_combined = similar(grid_densities[1], Int)
    grid_singles = [similar(grid_combined) for _ = 1:n_groups]

    sites = [
        Observable(
            Set(
                CartesianIndices(size(grid_combined))[rand(
                    1:length(grid_combined),
                    initial_sites,
                )],
            ),
        ) for _ = 1:n_groups
    ]

    return grid_combined, grid_singles, sites
end

function Makie.plot!(sm::StippleMap; kwargs...)

    grid_densities = sm.densities
    n_groups = length(to_value(grid_densities))

    grid_combined, grid_singles, sites =
        setup(to_value(grid_densities), to_value(sm.n_initial_sites))

    # init the grids
    jfa_voronoi_gpu!(grid_combined, collect(reduce(union, to_value.(sites))))
    [jfa_voronoi_gpu!(grid_singles[k], collect(to_value(sites[k]))) for k = 1:n_groups]


    on(sm.one_more_iteration; update = true) do trg
        _sites =
            Set.(linearized_sites.(collect.(to_value.(sites)), Ref(size(grid_combined))))
        @show trg
        for t = 1:trg
            one_iteration!(
                grid_combined,
                grid_singles,
                _sites,
                to_value(grid_densities);
                markersize = to_value(sm.markersize),
                hysteresis_a = to_value(sm.hysteresis),
                density_factor = to_value(sm.density_factor),
            )
        end
        for k = 1:n_groups
            sites[k][] =
                Set(linearized_sites.(collect(_sites[k]), Ref(size(grid_combined))))
        end
        #notify.(sites)

    end




    for k = 1:n_groups
        order = @lift(.-Matrix($(grid_densities)[k])[collect($(sites[k]))])
        site_obs = @lift(convert.(Point3, (collect($(sites[k]))), $order))
        #@debug typeof(site_obs) typeof(to_value(site_obs[1]))
        c = sm.colors
        kw = try
            (; color = @lift($c[k]))
        catch
            (;)
        end
        scatter!(
            sm,
            site_obs;
            markerspace = :data,
            markersize = sm.markersize,
            fxaa = true,
            kw...,
        )

    end
    return sm
end

Base.convert(Point2, x::CartesianIndex) = Point2(x[1], x[2])
Base.convert(Point3, x::CartesianIndex, y) = Point3(x[1], x[2], y)
