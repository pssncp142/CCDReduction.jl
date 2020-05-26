
import Base.Broadcast: broadcastable
import Base
using FITSIO

#=============================================================================#
# CCD Lazy
struct CCDLazy
    f_name::String
    length::Int
    size::Tuple{Int,Int}
    function CCDLazy(f_name::String)
        hdu = FITS(f_name)
        len = length(hdu[1])
        sz = size(hdu[1])
        return new(f_name, len, sz)
    end
end

broadcastable(im::CCDLazy) = broadcastable(convert(CCDImage, im))

Base.:length(im::CCDLazy) = im.length
Base.:size(im::CCDLazy) = im.size

# end CCDLazy
#=============================================================================#
# CCDImage Definitions

mutable struct CCDImage <: AbstractArray{Pixel,2}
    data::Array{Float64,2}
    error::Array{Float64,2}
    mask::Array{Int8,2}
    function CCDImage(data, error, mask=-1)
        shape = size(data)
        if mask == -1
            mask = zeros(shape)
        end
        new(data, error, mask)
    end
end

CCDImage(ndx1::Int64, ndx2::Int64) = CCDImage(zeros(ndx1, ndx2), zeros(ndx1,ndx2), zeros(ndx1,ndx2))

CCDImage(f_name::String) = CCDLazy(f_name)

function Base.:convert(T::Type{CCDImage}, im::CCDLazy)
    hdu = FITS(im.f_name)
    data = read(hdu[1])
    shape = size(hdu[1])
    return CCDImage(data, ones(shape), zeros(shape))
end

# end CCDImage
#=============================================================================#
# Array Functionality

Base.:BroadcastStyle(::Type{<:CCDImage}) = Broadcast.ArrayStyle{CCDImage}()
Base.:eltype(::Type{CCDImage}) = Pixel

function Base.:similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{CCDImage}}, ::Type{ElType}) where ElType
    return CCDImage(zeros(axes(bc)), zeros(axes(bc)), zeros(axes(bc)))
end

function Base.:getindex(im::CCDImage, I::Vararg{UnitRange{Int64}})
    res = CCDImage(I[1][end]-I[1][1]+1, I[2][end]-I[2][1]+1)
    for i in I[1]
        for j in I[2]
            res[i-I[1][1]+1,j-I[2][1]+1] = im[i,j]
        end
    end
    return res
end

Base.:getindex(im::CCDImage, I::Vararg{Int64,2}) = Pixel(im.data[I...], im.error[I...], im.mask[I...])

function Base.:setindex!(im::CCDImage, pix::Pixel, I::Vararg{Int64,2})
    if nfields(I) == 1
        im.data[I[1]] = pix.data
        im.error[I[1]] = pix.error
        im.mask[I[1]] = pix.mask
    else
        im.data[I[1],I[2]] = pix.data
        im.error[I[1],I[2]] = pix.error
        im.mask[I[1],I[2]] = pix.mask
    end
end

Base.:setindex!(im::CCDImage, val::Number, I::Vararg{Int64,2}) = setindex!(im, Pixel(val), I...)

function Base.:iterate(im::CCDImage, ndx=1)
	return ndx > length(im) ? nothing : (Pixel(im.data[ndx], im.error[ndx], im.mask[ndx]), ndx+1)
end

Base.:length(im::CCDImage) = length(im.data)
Base.:size(im::CCDImage) = Base.size(im.data)

# end Array Functionality
#=============================================================================#
