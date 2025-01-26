
"""
translates sites from (x,y) format to indices and the other way round

# todo make consistent
"""
linearized_sites(p::Int, sz_grid) = CartesianIndices(sz_grid)[p]
function linearized_sites(sites, sz_grid)
    function _mycheck(x)
        if x[1] <= 0 || x[2] <= 0
            return 0
        end
        LinearIndices(sz_grid)[x[1], x[2]]
    end
    map(x -> _mycheck(x), sites)#
end
#(getindex.(sites, 2) .- 1) * sz_grid[1] .+ getindex.(sites, 1)
# LinearIndices(sz_grid)



function run_iterations!(
    grid_combined,
    sites,
    sites_grouping,
    density_interpolators;
    n_iter = 10,
    kwargs...,
)
    n_conditions = length(density_interpolators)

    # initialize individual voronoi grids
    grid_single = Vector{Any}(undef, n_conditions)
    for k = 1:n_conditions
        grid_single[k] = similar(grid_combined)
        grid_single[k] .= 0
    end

    # generate the sparse active-site-vector to be modified
    p_all = linearized_sites(sites, size(grid_combined))
    sites_sets = [Set(p_all[sites_grouping.==k]) for k = 1:n_conditions]


    jfa_voronoi_gpu!(grid_combined, collect(reduce(union, sites_sets)))
    [jfa_voronoi_gpu!(grid_single[k], collect(sites_sets[k])) for k = 1:n_conditions]

    for it = 1:n_iter # iteration
        # init sparse vector indicating which point is active
        @show "iteration $it, n-points $(sum(length.(sites_sets)))"
        one_iteration!(
            grid_combined,
            grid_single,
            sites_sets,
            density_interpolators;
            kwargs...,
        )
        # calculate areas
    end
    return grid_single, sites_sets
end



function one_iteration!(
    grid_combined,
    grid_single,
    sites_sets,
    density_interpolators;
    kwargs...,
)
    n_conditions = length(density_interpolators)


    split_by_density!(
        sites_sets,
        grid_combined,
        grid_single,
        density_interpolators;
        kwargs...,
    )

    jfa_voronoi_gpu!(grid_combined, collect(reduce(union, sites_sets)))
    [jfa_voronoi_gpu!(grid_single[k], collect(sites_sets[k])) for k = 1:n_conditions]
end

function split_by_density!(
    sites_sets,
    grid_combined,
    grid_single,
    grid_densities;
    threshold_low = 10e-7,
    threshold_high = 0.0000005,
)
    n_conditions = length(grid_single)
    sites_combined = collect(reduce(union, sites_sets))
    areas_combined, centroids_combined = grid_features(grid_combined, sites_combined)

    centroids_single = Array{Any}(undef, n_conditions)
    centroids_single_linearized = Array{Any}(undef, n_conditions)
    areas_single = Array{Any}(undef, n_conditions)
    for (k, grid) in enumerate(grid_single)
        p_nearest = grid[linearized_sites.(sites_combined, Ref(size(grid)))]
        areas_single[k], centroids_single[k] = grid_features(grid, p_nearest)##grid_features(grid,p_nearest)

        # need this helper to do a double broadcast
        function _my_round(x)
            Int.(round.(x))
        end
        # round and linearize
        centroids_single_linearized[k] =
            linearized_sites(_my_round.(centroids_single[k]), size(grid_combined))
    end
    remove = 0
    split = 0
    for (ix, p) in enumerate(sites_combined)
        #---
        # calculate the required centroids / areas
        #centroid_combined, area_combined = grid_features(grid_combined,p)
        #calculate_density(density_interpolator, centroids_single[k][ix])
        #@show ix, centroids_single[1][ix]

        densities = Vector{Float64}(undef, n_conditions)
        for k = 1:n_conditions
            centroid_index = centroids_single_linearized[k][ix]
            if centroid_index <= 0 || centroid_index > *(size(grid_combined)...)
                densities[k] = 0.0
                # @warn "out of bound centroid $(centroids_single[k][ix]) found in condition=$k, ix=$ix"
                continue
            end
            densities[k] = grid_densities[k][centroid_index] * areas_single[k][ix]
        end

        #@show densities
        density_argmax = argmax(densities)
        density_selected = densities[density_argmax]
        centroid_selected = centroids_single[density_argmax][ix]

        probability = sum(densities) - density_selected
        #@debug probability
        centroid_new = (probability < rand()) ? centroid_selected : centroids_combined[ix]

        #@debug "centroids sel/comb" centroid_selected centroids_combined[ix]
        #@debug density_selected
        [delete!(sites_sets[k], p) for k = 1:length(densities)]
        #if centroid_new in sites)
        if density_selected < threshold_low
            remove += 1
            #@show "let's only remove..."
            # remove it
            #delete!(sites_sets[density_argmax],p)
            #[delete!(sites_sets[k],p) for k = 1:length(densities)]

        elseif density_selected > threshold_high
            #@show "splitting it"
            split += 1

            # split the cell in two
            Δmove = (rand(2) .- 0.5) .* size(grid_combined) * 0.01
            for _tmp in eachcol([Δmove .-Δmove])
                centroid_moved = move_centroid(centroid_new, _tmp, size(grid_combined))

                #LinearIndices(size(grid_combined))[centroid_new.+_tmp]#
                #p_new = linearized_sites([centroid_new.+_tmp],size(grid_combined))[1] # within 1% of the centroid_new
                p_new_int =
                    LinearIndices(size(grid_combined))[centroid_moved[1], centroid_moved[2]]

                #@debug "split centroids:" centroid_new,linearized_sites(p_new_int,size(grid_combined))
                #@debug "pushing " p_new_int
                push!(sites_sets[density_argmax], p_new_int)
                #@show centroid_new.+_tmp
                #@show p_new_int
                #@show linearized_sites(p_new_int,size(grid_combined))
            end

            #push!(sites_sets[density_argmax],Int(round(linearized_sites([centroid_new],size(grid_combined))[1])))

        else
            # keep it around, move to selected
            centroid_new = Int.(round.(centroid_new))
            centroid_new = min.(max.(1, centroid_new), size(grid_combined))
            push!(
                sites_sets[density_argmax],
                LinearIndices(size(grid_combined))[centroid_new[1], centroid_new[2]],
            )
            #push!(sites_sets[density_argmax],Int(round(linearized_sites([centroid_new],size(grid_combined))[1])))

        end


    end # end p-loop
    @debug "remove | split | total" remove, split, length(sites_combined)
end


function move_centroid(centroid_new, _tmp, sz_grid)
    centroid_moved = Int.(round.(centroid_new .+ _tmp))
    centroid_moved = min.(max.(1, centroid_moved), sz_grid)
end
