export
    AbstractParticle,
    AbstractMParticle,
    Particle,
    MParticle,
    MVParticle,
    Atom

abstract type AbstractParticle{D,T} end

abstract type AbstractMParticle{D,T} <: AbstractParticle{D,T} end

function Base.show(io::IO, p::T) where T<:AbstractParticle
    print(io, "$(T) $(p.type) (r = $(p.r))")
end # function

struct Particle{D,T} <: AbstractParticle{D,T}
    type::String #identifier
    r::SVector{D,T} #position
end # struct

Particle(; type = "", r = zeros(SVector{3})) = Particle(type, r)
Particle(r::SVector) = Particle(r = r)
Particle(r) = Particle(r = SVector{length(r)}(r))

struct Atom <: AbstractParticle{3,Float64}
    type::String
    resid::Int
    resname::String
    m::Float64
    q::Float64
    σ::Float64
    ϵ::Float64
    r::SVector{3,Float64}
    v::SVector{3,Float64}
end # struct

function Atom(;
              type = "",
              resid = 0,
              resname = "",
              m = 0.0,
              q = 0.0,
              σ = 0.0,
              ϵ = 0.0,
              r = zeros(SVector{3}),
              v = zero(r))
    return Atom(type, resid, resname, m, q, σ, ϵ, r, v)
end # function
    
function Base.show(io::IO, a::Atom)
    print(io, "Atom $(a.type) (r = $(a.r); v = $(a.v))")
end # function


struct MParticle{D,T} <: AbstractMParticle{D,T}
    type::String
    r::MVector{D,T}
end # struct

MParticle(; type = "", r = zeros(MVector{3})) = MParticle(type, r)
MParticle(r::MVector) = MParticle(r = r)
MParticle(r) = MParticle(r = MVector{length(r)}(r))


struct MVParticle{D,T} <: AbstractMParticle{D,T}
    type::String
    r::MVector{D,T}
    v::MVector{D,T}
end # struct

MVParticle(; type = "", r = zeros(MVector{3}), v = zeros(MVector{3})) = MVParticle(type, r, v)
MVParticle(r::MVector) = MVParticle(r = r)
MVParticle(r) = MVParticle(r = MVector{length(r)}(r))
MVParticle(r::MVector, v::MVector) = MVParticle(r = r, v = v)
MVParticle(r, v) = MVParticle(r = MVector{length(r)}(r), v = MVector{length(v)}(v))


