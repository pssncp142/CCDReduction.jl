module CCDReduction

using Statistics

export subtract_bias,
       subtract_bias!,
       subtract_overscan,
       subtract_overscan!,
       CCDImage,
       Pixel

include("methods.jl")
include("pixel.jl")
include("CCDImage.jl")

end
