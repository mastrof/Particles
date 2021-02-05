export
    readgro


@doc raw"""
    groline(s::AbstractString, type = :particle)

Parse a single line from a .gro file and return a `Particle` or `Atom` object as defined by the `type` argument.
"""
function groline(s::AbstractString, type = :particle)
    resid = parse(Int, SubString(s, 1:5))
    resname::String = SubString(s, 6:10) |> strip
    atomname::String = SubString(s, 11:15) |> strip
    atomid = parse(Int, SubString(s, 16:20))
    r = parse.(Float64, SubString(s, 21:44) |> split)
    if type == :particle
        return Particle(type = atomname, r = SVector{3}(r))
    elseif type == :atom
        return Atom(type = atomname, resid = resid,
                    resname = resname, r = r)
    end # if 
end # function


@doc raw"""
    readgro(filename::AbstractString)

Read a .gro file `filename` and return the number of atoms, simulation box size(only cubic boxes supported) and atomic configuration.
"""
function readgro(filename::AbstractString)
    title = readlines(filename, 1)[1]
    natoms = parse(Int, readlines(filename, 1, 1)[1])
    config = [groline(s) for s in readlines(filename, natoms, 2)]
    boxstring = split(readlines(filename, 1, natoms+2)[1], r" +")
    box = [parse(Float64, l) for l in filter(!isempty, boxstring)]
    return natoms, SVector{3}(box), config
end # function
