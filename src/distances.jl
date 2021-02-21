export
    distancevector,
    distance,
    msd,
    rmsd

function distancevector(p1::AbstractVector, p2::AbstractVector)
    return p2 .- p1
end # function

function distancevector(p1::AbstractVector, p2::AbstractVector, box)
    δr = (p2 .- p1) ./ box
    δr = (δr .- round.(δr)) .* box
    return δr
end # function

function distancevector(p1::AbstractParticle, p2::AbstractParticle)
    return p2.r .- p1.r
end # function

function distancevector(p1::AbstractParticle, p2::AbstractParticle, box)
    δr = (p2.r .- p1.r) ./ box
    δr = (δr .- round.(δr)) .* box
    return δr
end # function

distance(p1, p2) = norm(distancevector(p1, p2))
distance(p1, p2, box) = norm(distancevector(p1, p2, box))


function msd(trajectory::T) where {S<:Real,T<:AbstractArray{S,3}}
    D, nparticles, nsteps = size(trajectory)
    δ = zeros(nparticles)
    meansquaredisp = zeros(nsteps)
    for t in 2:nsteps
        for i in 1:nparticles
            u = @view trajectory[:, i, t]
            v = @view trajectory[:, i, t-1]
            δ[i] = distance(u, v)
        end # for
        meansquaredisp[t] = sum(abs2.(δ)) / nparticles
    end # for
    return cumsum(meansquaredisp)
end # function

function msd(trajectory::T) where {S<:AbstractParticle,T<:AbstractArray{S,2}}
    nparticles, nsteps = size(trajectory)
    δ = zeros(nparticles)
    meansquaredisp = zeros(nsteps)
    for t in 2:nsteps
        for i in 1:nparticles
            u = @view trajectory[i, t]
            v = @view trajectory[i, t-1]
            δ[i] = distance(u, v)            
        end # for
        meansquaredisp[t] = sum(abs2.(δ)) / nparticles
    end # for
    return cumsum(meansquaredisp)
end # function

function msd(trajectory::T) where {S<:Real,T<:AbstractArray{S,2}}
    D, nsteps = size(trajectory)
    δ = zeros(nsteps)
    for t in 2:nsteps
        u = @view trajectory[:, t]
        v = @view trajectory[:, t-1]
        δ[t] = distance(u, v)
    end # for
    squaredisp = cumsum(abs2.(δ))
    return squaredisp
end # function

function msd(trajectory::T) where {S<:AbstractParticle,T<:AbstractArray{S,1}}
    nsteps = size(trajectory, 1)
    δ = zeros(nsteps)
    for t in 2:nsteps
        u = @view trajectory[t]
        v = @view trajectory[t-1]
        δ[t] = distance(u, v)
    end # for
    squaredisp = cumsum(abs2.(δ))
    return squaredisp
end # function

rmsd(trajectory) = sqrt.(msd(trajectory))
