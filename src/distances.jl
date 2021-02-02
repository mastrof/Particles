export
    distancevector,
    distance

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
