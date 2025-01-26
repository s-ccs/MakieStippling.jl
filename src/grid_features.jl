

"""
returns area and centroids for all elements in sites_set.
"""
function grid_features(A::AbstractMatrix{<:Integer}, sites_set)
    device = get_backend(A)
    # Define the number of bins for the histogram
    num_bins = length(sites_set)

    # Create arrays to store the histogram counts
    hist_counts = similar(A, Int32, num_bins)
    hist_counts .= 0

    # Create arrays to store the sum of x-coordinates for each bin
    hist_x_sum = similar(A, Int32, num_bins)
    hist_x_sum .= 0
    # Create arrays to store the sum of y-coordinates for each bin
    hist_y_sum = similar(A, Int32, num_bins)
    hist_y_sum .= 0

    sites_array = similar(A, Int32, num_bins)
    sites_array = typeof(sites_array)(Int32.(collect(sites_set)))
    # Define the kernel
    @kernel function histogram_kernel!(A, hist_counts, hist_x_sum, hist_y_sum, sites_set)
        i = @index(Global, Linear)
        #x = mod(i - 1, size(A, 2)) + 1
        #y = cld(i, size(A, 2))
        sz_grid = size(A)
        coords = CartesianIndices(sz_grid)

        x = coords[i][1]
        y = coords[i][2]

        bin_idx = findfirst(==(A[x, y]), sites_set)
        if !isnothing(bin_idx)
            Atomix.@atomic hist_counts[bin_idx] += 1
            KernelAbstractions.@atomic hist_x_sum[bin_idx] += Int32(x)
            KernelAbstractions.@atomic hist_y_sum[bin_idx] += Int32(y)
            #=
            =#
        end
    end
    # @show "launch hist"
    # Launch the kernel
    event = histogram_kernel!(device, 128)(
        A,
        hist_counts,
        hist_x_sum,
        hist_y_sum,
        sites_array;
        ndrange = size(A),
    )
    KernelAbstractions.synchronize(device)


    # Calculate the centroids
    centroids_x = similar(A, Float32, num_bins)
    centroids_y = similar(A, Float32, num_bins)

    @kernel function centroid_kernel!(
        hist_counts,
        hist_x_sum,
        hist_y_sum,
        centroids_x,
        centroids_y,
    )
        i = @index(Global, Linear)
        if hist_counts[i] > 0
            centroids_x[i] = hist_x_sum[i] / hist_counts[i]
            centroids_y[i] = hist_y_sum[i] / hist_counts[i]
        end
    end
    #    @show "launch centroid"
    # Launch the centroid kernel
    event = centroid_kernel!(device, 128)(
        hist_counts,
        hist_x_sum,
        hist_y_sum,
        centroids_x,
        centroids_y;
        ndrange = (num_bins,),
    )
    KernelAbstractions.synchronize(device)

    return Vector(hist_counts) ./ *(size(A)...),
    vec2tuple.(Vector(centroids_x), Vector(centroids_y))
end
vec2tuple(x, y) = (x, y)



#=
naive CPU version

function grid_features(grid,id::Int)
    ix = findall(x->x==id,grid)
    centroid = Tuple(sum(ix))./length(ix)
    area = length(ix)
return centroid,area
end
=#
