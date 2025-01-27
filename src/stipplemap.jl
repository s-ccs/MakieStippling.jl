@recipe(StippleMap, densities) do scene
    Attributes(markersize = 10, hysteresis = 0.6)
end

function MakieCore.plot!(sm::StippleMap)
    grid_densities = sm.densities

    scatter!.(sm, sites, markerspace = :data, markersize = sm.markersize)
    return sm
end
