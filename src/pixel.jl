import Base

#=============================================================================#
# Pixel Struct

struct Pixel <: Real
    data::Float64
    error::Float64
    mask::Int8
end

# Init from Scalars
Pixel(val::Number) = Pixel(val, 0, 0)

# val+/-error(*)
Base.:repr(pix::Pixel) = string(pix.data,"+/-",pix.error, convert(Bool, pix.mask) ? "*" : "")
Base.:show(io::IO, pix::Pixel) = print(io, repr(pix))

#=============================================================================#
# Pixel Arithmetic

function Base.:-(pix::Pixel)
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
    error = abs(data)*sqrt((pix1.error/pix1.data)^2 + (pix2.error/pix2.data)^2)
    mask = pix1.mask | pix2.mask
    return Pixel(data, error, mask)
end

function divide(pix1::Pixel, pix2::Pixel)
    data = pix1.data / pix2.data
    error = abs(data)*sqrt((pix1.error/pix1.data)^2 + (pix2.error/pix2.data)^2)
    mask = pix1.mask | pix2.mask
    return Pixel(data, error, mask)
end

function Base.:sqrt(pix::Pixel)
    data = sqrt(pix.data)
    error = abs(data*0.5*pix.error/pix.data)
    return Pixel(data, error, pix.mask)
end

# Assumes val has no uncertainty
function Base.:^(pix::Pixel, val::Pixel)
    data = pix.data^val.data
    error = abs(data*val.data*pix.error/pix.data)
    return Pixel(data, error, pix.mask)
end

Base.:+(a::Pixel, b::Pixel) = add(a,b)
Base.:-(a::Pixel, b::Pixel) = subtract(a,b)
Base.:*(a::Pixel, b::Pixel) = multiply(a,b)
Base.:/(a::Pixel, b::Pixel) = divide(a,b)

Base.:>(a::Pixel, b::Pixel) = a.data > b.data ? true : false
Base.:<(a::Pixel, b::Pixel) = a.data < b.data ? true : false

# This enables all operations but totally ignores uncertainty as a result.
#Base.:AbstractFloat(pix::Pixel) = pix.data

Base.:convert(::Type{Pixel}, x::T) where T <: Union{Integer,AbstractFloat} = Pixel(x)

Base.:promote_rule(::Type{Pixel}, ::Type{T}) where T <: Union{Integer,AbstractFloat} = Pixel

Base.:one(::Type{Pixel}) = Pixel(1)
Base.:zero(::Type{Pixel}) = Pixel(0)

# end Pixel Arithmetic
#=============================================================================#

