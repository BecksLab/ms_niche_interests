using Graphs

"""
randommodel(species::Int64, L::Int64)

    Return a network of randomly assembled interactions according to
    the Erdős-Rényi model.

    This is essentially a wrapper for the `erdos_renyi` function from
    Graphs.jl and just packages it into an adj. matrix.

    species = number of species
    L = number of links

    #### References

    Erdős, Paul, and Alfréd Rényi. 1959. “On Random Graphs I.” 
    Publicationes Mathematicae. https://doi.org/10.5486/PMD.1959.6.3-4.12.

    Graphs.jl TODO
"""
function random(species::Any, Co::Int64)

    L = floor(Int, Co * (length(species)^2))

    S = length(species)
    N = erdos_renyi(S, L)

    # empty matrix
    edges = zeros(Bool, (S, S))

    for i = 1:S
        if length(N.fadjlist[i]) > 0
            for j in eachindex(N.fadjlist[i])
                edges[i, j] = 1
            end
        end
    end

    return edges
end
