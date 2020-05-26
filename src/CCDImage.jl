
import Base: +,-,*,/,getindex,setindex!,length,size,iterate,IteratorSize,BroadcastStyle,similar,eltype,show,showarg,repr,convert,one,zero,promote_rule
import Base.Broadcast: broadcastable

import Statistics: sum,mean,median

using FITSIO

#=============================================================================#
# Pixel Struct

struct Pixel <: Number
    data::Float64
    error::Float64
    mask::Int8
end

# val+/-error(*)
repr(pix::Pixel) = string(repr(pix.data),"+/-",repr(pix.error), convert(Bool, pix.mask) ? "*" : "")

show(io::IO, pix::Pixel) = print(io, repr(pix))

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

length(im::CCDLazy) = im.length
size(im::CCDLazy) = im.size

# end Pixel Struct
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

function convert(T::Type{CCDImage}, im::CCDLazy)
    hdu = FITS(im.f_name)
    data = read(hdu[1])
    shape = size(hdu[1])
    return CCDImage(data, ones(shape), zeros(shape))
end

# end CCDImage
#=============================================================================#
# Pixel Arithmetic
# This is cleaner...

function -(pix::Pixel)
    data = -pix.data
    error = pix.error
    mask = pix.mask
    return Pixel(data, error, mask)
end

function add(pix1::Pixel, pix2::Pixel)
    data = pix1.data + pix2.data
    error = sqrt(pix1.error^2 + pix2.error^2)
    mask = pix1.mask | pix2.mask
    return Pixel(data, error, mask)
end

function subtract(pix1::Pixel, pix2::Pixel)
    data = pix1.data - pix2.data
    error = sqrt(pix1.error^2 + pix2.error^2)
    mask = pix1.mask | pix2.mask
    return Pixel(data, error, mask)
end

function multiply(pix1::Pixel, pix2::Pixel)
    data = pix1.data * pix2.data
    error = sqrt((pix2.data * pix1.error)^2 + (pix1.data * pix2.error)^2)
    mask = pix1.mask | pix2.mask
    return Pixel(data, error, mask)
end

function divide(pix1::Pixel, pix2::Pixel)
    data = pix1.data / pix2.data
    error = sqrt((pix2.data * pix1.error)^2 + (pix1.data * pix2.error)^2)
    mask = pix1.mask | pix2.mask
    return Pixel(data, error, mask)
end

(+)(a::Pixel, b::Pixel) = add(a,b)
(-)(a::Pixel, b::Pixel) = subtract(a,b)
(*)(a::Pixel, b::Pixel) = multiply(a,b)
(/)(a::Pixel, b::Pixel) = divide(a,b)

convert(::Type{Pixel}, x::Number) = Pixel(x, 0, 0)
convert(::Type{Pixel}, pix::Pixel) = pix

promote_rule(::Type{Pixel}, ::Type{<:Number}) = Pixel

one(::Type{Pixel}) = Pixel(1, 0, 0)
zero(::Type{Pixel}) = Pixel(0, 0, 0)

# end Pixel Arithmetic
#=============================================================================#
# Array Functionality

BroadcastStyle(::Type{<:CCDImage}) = Broadcast.ArrayStyle{CCDImage}()
eltype(::Type{CCDImage}) = Pixel

function similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{CCDImage}}, ::Type{ElType}) where ElType
    return CCDImage(ones(axes(bc)), ones(axes(bc)), ones(axes(bc)))
end

function getindex(im::CCDImage, I::Vararg{T,2}) where T Union{Int64,UnitRange{Int64}}
    rn1 = 1:1
    rn2 = 1:1
    if typeof(I[1]) == Int64 && typeof(I[2]) == Int64
        return Pixel(im.data[I[1],I[2]], im.error[I[1],I[2]], im.mask[I[1],I[2]])
    elseif typeof(I[1]) == Int64
        rn1 = I[1]:I[1]
    elseif typeof(I[2]) == Int64
        rn2 = I[2]:I[2]
    else
        rn1 = I[1]
        rn2 = I[2]
    end
    res = CCDImage(rn1[end]-rn1[1]+1, rn2[end]-rn2[1]+1)
    for i in rn1
        for j in rn2
            res[i-rn1[1]+1,j-rn2[1]+1] = im[i,j]
        end
    end
    return res
end

function setindex!(im::CCDImage, pix::Pixel, I::Vararg{T,2}) where T Union{Int64}
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

function iterate(im::CCDImage, ndx=1)
	return ndx > length(im) ? nothing : (Pixel(im.data[ndx], im.error[ndx], im.mask[ndx]), ndx+1)
end

function length(im::CCDImage)
    return length(im.data)
end

function size(im::CCDImage)
    return size(im.data)
end

# end Array Functionality
#=============================================================================#
# Simple Statistics (Just as an example)

function sum(im::CCDImage)
    return sum(im.data[im.mask .== false])
end

function mean(im::CCDImage)
    return mean(im.data[im.mask .== false])
end

function median(im::CCDImage)
    return median(im.data[im.mask .== false])
end

# end Simple Statistics
#=============================================================================#
