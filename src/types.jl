export
    AbstractParticle,
    AbstractMParticle,
    Particle,
    MParticle,
    MVParticle,
    Atom

@doc raw"""
    AbstractParticle{D,T}

General interface for particle-like data.
`D` is the dimensionality of the `AbstractParticle` (the number of coordinates that describe its position).
`T` is the type of the `AbstractParticle` space (e.g., `Int` for particles on a lattice, `Float64` for real space, or even `ComplexF64`).

`AbstractParticle` is on top of the hierarchy
    `AbstractParticle` <: `AbstractMParticle`
`AbstractMParticle` uses mutable structs, most convenient for particles whose position can change throughout a program.
"""
abstract type AbstractParticle{D,T} end

abstract type AbstractMParticle{D,T} <: AbstractParticle{D,T} end

function Base.show(io::IO, p::T) where T<:AbstractParticle
    print(io, "$(T) $(p.type) (r = $(p.r))")
end # function

Base.zero(p::T) where {T<:AbstractParticle} = T(type = p.type)



struct Particle{D,T} <: AbstractParticle{D,T}
    type::String #identifier
    r::SVector{D,T} #position
end # struct

Particle(type, r::AbstractVector) = Particle(type, SVector{length(r)}(r))
Particle{D,T}(; type = "", r = zeros(SVector{D,T})) where {D,T} =
    Particle(type, r)
Particle(; type = "", r = zeros(SVector{3})) = Particle(type, r)
Particle(r::SVector) = Particle(r = r)
Particle(r::AbstractVector) = Particle(r = SVector{length(r)}(r))


struct MParticle{D,T} <: AbstractMParticle{D,T}
    type::String #identifier
    r::MVector{D,T} #position
end # struct

MParticle(type, r::AbstractVector) =
    MParticle(type, MVector{length(r)}(r))
MParticle{D,T}(; type = "", r = zeros(MVector{D,T})) where {D,T} =
    MParticle(type, r)
MParticle(; type = "", r = zeros(MVector{3})) = MParticle(type, r)
MParticle(r::MVector) = MParticle(r = r)
MParticle(r::AbstractVector) = MParticle(r = MVector{length(r)}(r))


struct MVParticle{D,T} <: AbstractMParticle{D,T}
    type::String #identifier
    r::MVector{D,T} #position
    v::MVector{D,T} #velocity
end # struct

MVParticle(type, r::AbstractVector, v::AbstractVector) =
    MVParticle(type, MVector{length(r)}(r), MVector{length(v)}(v))
MVParticle{D,T}(; type = "", r = zeros(MVector{D,T}),
                v = zeros(MVector{D,T})) where {D,T} =
                    MVParticle(type, r, v)
#=
MVParticle{D}(; type = "", r = zeros(MVector{D,Float64}),
              v = zeros(MVector{D,Float64})) where D =
                  MVParticle{D,Float64}(type, r, v)
MVParticle(; type = "") = MVParticle{3,Float64}(type=type)
=#
function MVParticle(; type = "", r = nothing, v = nothing)
    if isnothing(r) && isnothing(v)
        MVParticle{3,Float64}(type = type)
    elseif !isnothing(r) && isnothing(v)
        D = length(r)
        T = eltype(r)
        MVParticle{D,T}(type = type, r = r, v = zeros(MVector{D,T}))
    elseif isnothing(r) && !isnothing(v)
        D = length(v)
        T = eltype(v)
        MVParticle{D,T}(type = type, r = zeros(MVector{D,T}), v = v)
    else
        MVParticle(type, r, v)
    end # if
end # function
####
MVParticle(r::MVector) = MVParticle(r = r, v = zeros(MVector{length(r)}))
MVParticle(r::AbstractVector) = MVParticle(MVector{length(r)}(r))
MVParticle(r::MVector, v::MVector) = MVParticle(r = r, v = v)
MVParticle(r::AbstractVector, v::AbstractVector) =
    MVParticle(r = MVector{length(r)}(r), v = MVector{length(v)}(v))



struct Atom <: AbstractParticle{3,Float64}
    type::String
    resid::Int
    resname::String
    m::Float64
    q::Float64
    ??::Float64
    ??::Float64
    r::SVector{3,Float64}
    v::SVector{3,Float64}
end # struct

function Atom(;
              type = "",
              resid = 0,
              resname = "",
              m = 0.0,
              q = 0.0,
              ?? = 0.0,
              ?? = 0.0,
              r = zeros(SVector{3}),
              v = zero(r))
    return Atom(type, resid, resname, m, q, ??, ??, r, v)
end # function
    
function Base.show(io::IO, a::Atom)
    print(io, "Atom $(a.type) (r = $(a.r); v = $(a.v))")
end # function
