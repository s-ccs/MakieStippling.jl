
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
    grid_densities;
    n_iter = 10,
    kwargs...,
)
    grid_singles = [similar(grid_combined) for _ = 1:length(grid_densities)]
    grid_singles = setup(grid_combined, grid_densities)
    function setup(grid_combined, grid_densities)
        n_conditions = length(grid_densities)

        # initialize individual voronoi grids
        grid_singles = Vector{Any}(undef, n_conditions)
        for k = 1:n_conditions
            grid_singles[k] = similar(grid_combined)
            grid_singles[k] .= 0
        end

        # generate the sparse active-site-vector to be modified
        #p_all = linearized_sites(sites, size(grid_combined))
        #sites_sets = [Set(p_all[sites_grouping.==k]) for k = 1:n_conditions]
        return grid_singles
    end

    jfa_voronoi_gpu!(grid_combined, collect(reduce(union, sites_sets)))
    [jfa_voronoi_gpu!(grid_singles[k], collect(sites_sets[k])) for k = 1:n_conditions]

    for it = 1:n_iter # iteration
        # init sparse vector indicating which point is active
        @show "iteration $it, n-points $(sum(length.(sites_sets)))"
        one_iteration!(grid_combined, grid_singles, sites_sets, grid_densities; kwargs...)
        # calculate areas
    end
    return grid_single, sites_sets
end



function one_iteration!(grid_combined, grid_single, sites_sets, grid_densities; kwargs...)
    n_conditions = length(grid_densities)


    split_by_density!(sites_sets, grid_combined, grid_single, grid_densities; kwargs...)

    jfa_voronoi_gpu!(grid_combined, collect(reduce(union, sites_sets)))
    [jfa_voronoi_gpu!(grid_single[k], collect(sites_sets[k])) for k = 1:n_conditions]
end

# need this helper to do a double broadcast
function _my_round(x)
    max.(1, Int.(round.(x)))
end

function split_by_density!(
    sites_sets,
    grid_combined,
    grid_single,
    grid_densities;
    threshold_low = 10e-7,
    threshold_high = 0.0000005,
    markersize = 5,
    hysteresis_a = 0.6,
    density_factor = 0.5,
)
    if !isnothing(hysteresis_a)
        markerarea = 2 * π * (markersize * 0.7)^2 * density_factor# *0.7 because a marker does not have a diameter of 1 by defaultm but rather ~0.7, diameter==1 => markersize * 1.4, but then we need radius again, thus we divide again by 2
        threshold_low = (1 - hysteresis_a) * markerarea
        threshold_high = (1 + hysteresis_a) * markerarea
        #        @debug threshold_low, threshold_high


    end
    n_conditions = length(grid_single)
    sites_combined = collect(reduce(union, sites_sets)) #x/y coordinates of all sites
    areas_combined, centroids_combined = grid_features(grid_combined, sites_combined) # voronoi features of all sites
    #centroids_combined_linearized =
    #    linearized_sites(_my_round.(centroids_combined), size(grid_combined))
    centroids_single = Array{Any}(undef, n_conditions)
    centroids_single_linearized = Array{Any}(undef, n_conditions)
    areas_single = Array{Any}(undef, n_conditions)
    lookup_dicts = Array{Any}(undef, n_conditions)

    _densities = Vector{Any}(undef, n_conditions)
    probabilities = Vector{Any}(undef, n_conditions)

    for (k, grid) in enumerate(grid_single) # go through each subvoronoi

        #p_nearest = Vector(grid[sites_combined]) # get the nearest centroid in the subvoronoi for all sites
        #lookup_closest[k] = Dict(sites_combined .=> p_nearest) # make a lookup table
        #un_p_nearest = (unique(p_nearest)) # multiple combined-points will map to

        #p_nearest_indices = findall(in(un_p_nearest), p_nearest) # find indices in un_p_nearest of each element in p_nearest

        #@debug un_p_nearest
        #@debug lookup_closest[k]
        areas_single[k], centroids_single[k] = grid_features(grid, sites_sets[k])
        lookup_dicts[k] = Dict(sites_sets[k] .=> 1:length(sites_sets[k]))
        #        _areas, _centroids = grid_features(grid, un_p_nearest)
        #ix = map(x -> lookup_dict[x], p_nearest) #lookup_dict[un_p_nearest]
        #areas_single[k] = _areas[ix]
        #centroids_single[k] = _centroids[ix]


        # round and linearize

        centroids_single_linearized[k] =
            linearized_sites(_my_round.(centroids_single[k]), size(grid_combined))#=Dict(
        un_p_nearest .=>
            linearized_sites(_my_round.(_centroids), size(grid_combined)),
        )=#
        _densities[k] =
            Vector(grid_densities[k][centroids_single_linearized[k]] .* areas_single[k]) # ./
        D = sum(
            reduce(
                hcat,
                [
                    grid_densities[k_o][centroids_single_linearized[k]] .* areas_single[k] for k_o = 1:length(grid_single) if k != k_o
                ],
            ),
            dims = 2,
        )[
            :,
            1,
        ]

        probabilities[k] = Vector(D)# .- _densities[k] # we excluded K in the sum before :)
    end


    # now find out the densities for the other layers, dependent on the current layer
    #=
    probabilities = Vector{Any}(undef, n_conditions)
    for k = 1:length(grid_single)
        _dens = reduce(
            hcat,
            [
                _densities[k_other][lookup_dicts[k_other][grid_other[sites_sets[k]]]]
                for (k_other, grid_other) in enumerate(grid_single) if k != k_other
            ],
        )
        probabilities[k] = sum(_dens, dims = 2)[:, 1] # because we never added the d_l to D, we dont need to subtract
    end
    =#
    # calculate densities for all
    #densities = reduce(hcat, [Vector(d) for d in _densities])
    remove = 0
    split = 0
    #densities_argmax = argmax(densities, dims = 2)[:, 1]
    #densities_argmax = argmin([
    #    findmin(euclidean.(centroids_single[k],centroids_combined)) for k = 1:length(grid_single)
    #])

    #@debug size(densities_argmax) size(densities)
    for k = 1:length(grid_single)
        for (p_ix, p) in enumerate(collect(sites_sets[k]))
            #density_argmax = densities_argmax[ix][2]

            #p_ix = lookup_dicts[k][p]
            # density_selected = densities[ix, density_argmax]
            #p_k = lookup_closest[density_argmax][p]

            #centroid_selected = centroids_single[density_argmax][ix]

            #probability = densities_sum[ix] - density_selected
            #@debug probability centroid_selected centroids_combined[ix]
            c_ix = findfirst(==(p), sites_combined)
            centroid_new =
                (probabilities[k][p_ix] > rand()) ? centroids_single[k][p_ix] :
                centroids_combined[c_ix]

            # probably unnecessary, a point shouldnt be in multiple sets
            [delete!(sites_sets[k], p) for k = 1:length(grid_single)]
            #if centroid_new in sites)
            if _densities[k][p_ix] < threshold_low
                remove += 1

            elseif _densities[k][p_ix] > threshold_high
                #@show "splitting it"
                split += 1

                # split the cell in two
                Δmove = (rand(2) .- 0.5) .* markerarea#markerarea * 10#size(grid_combined) #0.01

                for _tmp in eachcol([Δmove .-Δmove])
                    centroid_moved = move_centroid(centroid_new, _tmp, size(grid_combined))
                    p_new_int = LinearIndices(size(grid_combined))[
                        centroid_moved[1],
                        centroid_moved[2],
                    ]

                    push!(sites_sets[k], p_new_int)
                end

            else
                # keep it around, move to selected
                centroid_new = Int.(round.(centroid_new))
                centroid_new = min.(max.(1, centroid_new), size(grid_combined))
                push!(
                    sites_sets[k],
                    LinearIndices(size(grid_combined))[centroid_new[1], centroid_new[2]],
                )


            end


        end # end p-loop
    end
    @debug "remove | split | total" remove, split, length(sites_combined)
end


function move_centroid(centroid_new, _tmp, sz_grid)
    centroid_moved = Int.(round.(centroid_new .+ _tmp))
    centroid_moved = min.(max.(1, centroid_moved), sz_grid)
end
