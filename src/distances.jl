export
    distancevector,
    distance,
    unfold,
    msd,
    rmsd

@doc raw"""
    distancevector(p1::AbstractParticle, p2::AbstractParticle)
    distancevector(p1::AbstractVector, p2::AbstractVector)

Evaluates the vector distance of AbstractParticle `p2` from AbstractParticle `p1`.

---
    distancevector(p1::AbstractParticle, p2::AbstractParticle, box)
    distancevector(p1::AbstractVector, p2::AbstractVector, box)

Evaluates the vector distance of AbstractParticle `p2` from AbstractParticle `p1` with periodic boundary conditions in a box of size `box`.
`box` can be either a Float64 (cubic box) or an AbstractVector (one value for each dimension of the particle space).
_Currently only works with orthorhombic boxes._
"""
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

@doc raw"""
    distance(p1, p2)
    distance(p1, p2, box)

Evaluates the (Euclidean) distance of AbstractParticle p2 from p1.
See also `distancevector`.
"""
distance(p1, p2) = norm(distancevector(p1, p2))
distance(p1, p2, box) = norm(distancevector(p1, p2, box))



@doc raw"""
    msd(trajectory)
    msd(trajectory, box)

Evaluate the mean-square displacement (MSD) of an AbstractParticle or a collection thereof:
    MSD(t) = <|x(t) - x(0)|²> = Σᵢ|xᵢ(t) - xᵢ(0)|² / N
If `box` is specified, the MSD is evaluated with periodic boundary conditions (_only orthorhombic boxes_).
If periodic boundary conditions are used, the trajectory is unfolded under the assumption that no particle has crossed the box more than once between successive timesteps.

`trajectory` can be of the following types:
- T where {S<:Real,T<:AbstractArray{S,3}} => array of size (D,np,ns) with D = dimensionality, np = num particles, ns = num steps
- T where {S<:AbstractParticle,T<:AbstractArray{S,2}} => matrix of size (np, ns) with np = num particles, ns = num steps
- T where {S<:Real,T<:AbstractArray{S,2}} => matrix of size (D, ns) for single-particle trajectory
- T where {S<:AbstractParticle,T<:AbstractArray{S,1}} => vector of size ns for single-particle trajectory

All methods return a Vector{Float64} of length equal to `size(trajectory)[end]`.
"""
function msd(trajectory::T) where {S<:Real,T<:AbstractArray{S,3}}
    D, nparticles, nsteps = size(trajectory)    
    meansquaredisp = zeros(nsteps)
    for dt in 1:nsteps-1
        for t0 in 1:nsteps-dt
            for i in 1:nparticles
                u = @view trajectory[:, i, t0]
                v = @view trajectory[:, i, t0+dt]
                δ = distance(u, v)
                meansquaredisp[1+dt] += δ*δ
            end # for
        end # for
        meansquaredisp[1+dt] /= (nparticles * (nsteps-dt))
    end # for
    return meansquaredisp
end # function

function msd(trajectory::T) where {S<:AbstractParticle,T<:AbstractArray{S,2}}
    nparticles, nsteps = size(trajectory)    
    meansquaredisp = zeros(nsteps)
    for dt in 1:nsteps-1
        for t0 in 1:nsteps-dt
            for i in 1:nparticles
                u = trajectory[i, t0]
                v = trajectory[i, t0+dt]
                δ = distance(u, v)
                meansquaredisp[1+dt] += δ*δ
            end # for            
        end # for
        meansquaredisp[1+dt] /= (nparticles * (nsteps-dt))
    end # for
    return meansquaredisp
end # function

function msd(trajectory::T) where {S<:Real,T<:AbstractArray{S,2}}
    D, nsteps = size(trajectory)
    meansquaredisp = zeros(nsteps)
    for dt in 1:nsteps-1
        for t0 in 1:nsteps-dt
            u = @view trajectory[:, t0]
            v = @view trajectory[:, t0+dt]
            δ = distance(u, v)
            meansquaredisp[1+dt] += δ*δ
        end # for
        meansquaredisp[1+dt] /= (nsteps-dt)
    end # for    
    return meansquaredisp
end # function

function msd(trajectory::T) where {S<:AbstractParticle,T<:AbstractArray{S,1}}
    nsteps = size(trajectory, 1)
    meansquaredisp = zeros(nsteps)
    for dt in 1:nsteps-1
        for t0 in 1:nsteps-dt
            u = trajectory[t0]
            v = trajectory[t0+dt]
            δ = distance(u, v)
            meansquaredisp[1+dt] += δ*δ
        end # for
        meansquaredisp[1+dt] /= (nsteps-dt)
    end # for    
    return meansquaredisp
end # function

### msd with unfolding

function unfold(conf1::T, conf0::T, box) where {S<:AbstractParticle,T<:AbstractVector{S}}
    getdim(p::AbstractParticle{D,T}) where {D,T} = D
    D = getdim(first(conf1))
    nparticles = size(conf1, 1)
    unfolded = zero.(conf1)
    boxtype = typeof(box)
    for i in 1:nparticles
        ptype = conf1[i].type
        x = distancevector(zero(conf1[i]), conf1[i])
        dx = distancevector(conf0[i], conf1[i])
        newx = zeros(D)
        for μ in 1:D            
            L = boxtype <: AbstractArray ? box[μ] : box            
            sdx = sign(dx[μ])
            α = round(abs(dx[μ] / L))
            newx[μ] = abs(dx[μ]) > L/2 ? x[μ]-L*α*sdx : x[μ]
        end # for
        unfolded[i] = Particle(ptype, SVector{D}(newx))
    end # for
    return unfolded
end # function

function unfold(trajectory::T, box) where {S<:AbstractParticle,T<:AbstractMatrix{S}}
    nparticles, nsteps = size(trajectory)
    unfolded = zero.(trajectory)
    unfolded[:,1] .= trajectory[:,1]
    for t in 2:nsteps
        oldconf = unfolded[:,t-1]
        newconf = trajectory[:,t]
        unfolded[:,t] .= unfold(newconf, oldconf, box)
    end # for
    return unfolded
end # function

function unfold(trajectory::T, box) where {S<:AbstractParticle,T<:AbstractVector{S}}
    nsteps = size(trajectory, 1)
    unfolded = permutedims(zero.(trajectory))
    unfolded[:,1] = [trajectory[1]]
    for t in 2:nsteps
        oldconf = unfolded[:,t-1]
        newconf = [trajectory[t]]
        unfolded[:,t] = unfold(newconf, oldconf, box)
    end # for
    return unfolded[1,:]
end # function

function unfold(trajectory::T, box) where {S<:Real,T<:AbstractArray{S,3}}
    D, np, ns = size(trajectory)
    newtraj = [Particle(trajectory[:,i,t]) for i in 1:np, t in 1:ns]
    unfolded = unfold(newtraj,box)
    [unfolded[i,t].r[μ] for μ in 1:D, i in 1:np, t in 1:ns]
end # function

function unfold(trajectory::T, box) where {S<:Real,T<:AbstractArray{S,2}}
    D, ns = size(trajectory)
    newtraj = [Particle(trajectory[:,t]) for t in 1:ns]
    unfolded = unfold(newtraj, box)
    [unfolded[t].r[μ] for μ in 1:D, t in 1:ns]
end # function

function msd(trajectory, box)
    unfolded = unfold(trajectory, box)
    msd(unfolded)
end # function


rmsd(trajectory) = sqrt.(msd(trajectory))
rmsd(trajectory, box) = sqrt.(msd(trajectory, box))
