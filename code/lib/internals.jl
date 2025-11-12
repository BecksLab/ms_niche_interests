# dependencies
using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics


"""
  build_network(
    adj_mat::Matrix{Int64},
)

    Build an object of type `SpeciesInteractionNetwork{Unipartite{Symbol}, Binary{Bool}}`
    for an adjacency matrix of type `Matrix{Int64}`. Note that it is assumed that the species
    are not specified or retained and are specified as 'generic'. Additionally the matrix should
    be square and rows and columns are mapped.

    # Arguments
    - `adj_mat::Matrix{Int64}`: An adjacency matrix (S x S).

    # Returns
    - `SpeciesInteractionNetwork{Unipartite{Symbol}, Binary{Bool}}`: Network of type
    SpeciesInteractionNetwork{Unipartite{Symbol}, Binary{Bool}}.
    
"""
function build_network(adj_mat::Matrix{Int64})

    # checks
    if size(adj_mat)[1] != size(adj_mat)[2]
        throw(ArgumentError("Matrix is not square"))
    end

    # get spp richness (to create node list)
    spp_rich = size(adj_mat)[1]
    # create generic species list labelled 1:spp_rich
    spp_list = Symbol.(collect(1:1:spp_rich))

    # specify edges and nodes
    # make adj_mat Boolean
    edges = Binary(iszero.(adj_mat))
    nodes = Unipartite(spp_list)

    return SpeciesInteractionNetwork(nodes, edges)
    
end

"""
network_summary(N::SpeciesInteractionNetwork{<:Partiteness, <:Binary})

    returns the 'summary statistics' for a network
"""
function network_summary(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    _gen = SpeciesInteractionNetworks.generality(N)
    gen = collect(values(_gen))
    vul = collect(values(SpeciesInteractionNetworks.vulnerability(N)))
    ind_maxgen = findmax(gen)[2]

    L = links(N)
    S = richness(N)
    l_s = L / S

    tls = _trophic_level(N)

    D = Dict{Symbol,Any}(
        :richness => richness(N),
        :connectance => SpeciesInteractionNetworks.connectance(N),
        :complexity => complexity(N),
        :trophic_level => findmax(collect(values(tls)))[2],
        :generality => std(gen / l_s),
        :vulnerability => std(vul / l_s),
        # note for motif calcs one needs to exclude cannibalism (as per Stouffer)
        :S1 =>
            length(
                findmotif(motifs(Unipartite, 3)[1], remove_cannibals(N)),
            )/(richness(N)^2),
        :S2 =>
            length(
                findmotif(motifs(Unipartite, 3)[2], remove_cannibals(N)),
            )/(richness(N)^2),
        :S4 =>
            length(
                findmotif(motifs(Unipartite, 3)[4], remove_cannibals(N)),
            )/(richness(N)^2),
        :S5 =>
            length(
                findmotif(motifs(Unipartite, 3)[5], remove_cannibals(N)),
            )/(richness(N)^2),
    )

    return D
end

"""
remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Identifies and sets cannibalistic link to zero
"""
function remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    # get adj matrix
    S = richness(N)
    A = zeros(Bool, (S, S))
    for i in axes(A, 1)
        for j in axes(A, 2)
            if N.edges[i, j] == true
                A[i, j] = true
            end
        end
    end

    A[diagind(A)] .= false

    nodes = Unipartite(species(N))
    edges = Binary(A)
    network = SpeciesInteractionNetwork(nodes, edges)

    return network
end

"""
_diameter(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Calculates the diameter of a food web. Where diameter is the longest 
    shortest path between two nodes
"""
function _diameter(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    # extract species names
    spp = SpeciesInteractionNetworks.species(N)
    # empty vector for storing shortest path for each spp
    shortpath = zeros(Int64, length(spp))

    # get shortest path
    for i in eachindex(spp)
        shortpath[i] = length(shortestpath(N, spp[i]))
    end

    # return max shortest path
    return findmax(shortpath)[1]
end

"""
_trophic_level(N::SpeciesInteractionNetwork)

    Calculates the trophic level of all species in a network using the average 
    shortest path from the prey of species 𝑖 to a basal species

    Williams, Richard J., and Neo D. Martinez. 2004. “Limits to Trophic Levels 
    and Omnivory in Complex Food Webs: Theory and Data.” The American Naturalist 
    163 (3): 458–68. https://doi.org/10.1086/381964.
"""
function _trophic_level(N::SpeciesInteractionNetwork)

    sp = species(N)

    # dictionary for path lengths
    pls = Dict{Any,Any}()

    for i in eachindex(sp)

        # prey of spp i
        preys = collect(successors(N, sp[i]))

        # only continue if species has preys...
        if length(preys) > 0
            # for summing each path length
            pl_temp = 0
            for j in eachindex(preys)

                pl_temp += distancetobase(N, preys[j])

                pls[sp[i]] = 1 + (1/length(sp)) * pl_temp

            end
        else
            pls[sp[i]] = 1
        end
    end
    # return trophic level Dict
    return pls
end
