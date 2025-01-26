


struct DensityInterpolator
    mean_interpolator::Any
    std_interpolator::Any
    max_density::Any
end

function calculate_density(d::DensityInterpolator, point)
    #@show point
    mean = d.mean_interpolator(point[1])
    std = d.std_interpolator(point[1])
    pdf_value = pdf(Normal(mean, std), point[2])
    normalized_value = (pdf_value) / (d.max_density)
    return normalized_value
end
