
function jfa_voronoi_gpu!(grid, sites::AbstractVector{<:CartesianIndex}; kwargs...)
    sz_grid = size(grid)
    jfa_voronoi_gpu!(grid, linearized_sites(sites, sz_grid); kwargs...)
end

function jfa_voronoi_gpu!(
    grid,
    sites::Vector{Int};
    distance = euclidean,
    do_1plusjfa = false,
)
    grid .= 0
    grid[sites] = sites
    k = max(size(grid)...)

    if do_1plusjfa

        jfa_kernel(get_backend(grid), 128)(grid, 1, ndrange = size(grid))
        KernelAbstractions.synchronize(get_backend(grid))
    end

    while k > 1
        k = k ÷ 2 + k % 2
        jfa_kernel(get_backend(grid), 128)(grid, k, distance, ndrange = size(grid))
        KernelAbstractions.synchronize(get_backend(grid))
    end
end



@kernel function jfa_kernel(grid, k, distance)
    ix_grid = @index(Global) # could be slightly optimized with Global, Cartesian

    sz_grid = size(grid)
    coords = CartesianIndices(sz_grid)

    coords_grid = Tuple(coords[ix_grid])
    (x, y) = coords_grid

    for j in (-k, 0, k), i in (-k, 0, k)
        checkbounds(Bool, grid, x + i, y + j) || continue
        i == j == 0 && continue
        # get candidate at k-distanced grid location
        ix_siteq = grid[x+i, y+j]
        ix_siteq != 0 || continue
        # get current gridlocation value
        ix_sitep = grid[x, y]
        # if current location is not yet assigned, assign it with candidate
        if ix_sitep == 0
            grid[x, y] = ix_siteq
        elseif distance(Tuple(coords[ix_sitep]), coords_grid) >
               distance(Tuple(coords[ix_siteq]), coords_grid)
            # if candidate site q is closer to what is currently saved, use that
            grid[x, y] = ix_siteq
        end
    end


end
