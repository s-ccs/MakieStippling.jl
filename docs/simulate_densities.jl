function simulate_densities(; timepoints = 5, grid_sz, n_conditions = 2)

    times = collect(range(0, grid_sz[1], length = timepoints))
    density_interpolators = Vector{Any}(undef, n_conditions)
    max_density = 1
    for k = 1:n_conditions
        mean_data = Vector{Float64}(undef, timepoints)
        rng = MersenneTwister(k * 2 - 1)
        mean_data[1] = rand(rng)
        for t = 2:timepoints
            mean_data[t] = rand(rng) + mean_data[t-1] * 0.2
        end
        mean_data = mean_data .* grid_sz[1] / 4 .+ grid_sz[1] / 2

        std_data = abs.(rand(MersenneTwister(k * 2), timepoints)) .* grid_sz[2] * 0.1

        itp_mean = MakieStippling.Interpolations.interpolate(
            (times,),
            mean_data,
            MakieStippling.Interpolations.Gridded(MakieStippling.Interpolations.Linear()),
        )
        itp_std = MakieStippling.Interpolations.interpolate(
            (times,),
            std_data,
            MakieStippling.Interpolations.Gridded(MakieStippling.Interpolations.Linear()),
        )

        density_interpolators[k] =
            MakieStippling.DensityInterpolator(itp_mean, itp_std, max_density)

    end

    grid_densities = Vector(undef, n_conditions)
    for k = 1:n_conditions
        grid_densities[k] = Array{Float64}(undef, grid_sz)
        for x = 1:size(grid_densities[k], 1)
            for y = 1:size(grid_densities[k], 2)

                grid_densities[k][x, y] =
                    MakieStippling.calculate_density(density_interpolators[k], (x, y))
            end
        end
    end

    # normalize
    grid_densities = grid_densities ./ max(maximum.(grid_densities)...)
    return grid_densities
end
