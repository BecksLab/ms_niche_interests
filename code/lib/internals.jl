# dependencies
using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics


"""
  build_network(
    adj_mat::AbstractMatrix{Int},
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
function build_network(adj_mat::AbstractMatrix{Int})

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
    edges = Binary(.!iszero.(adj_mat))
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
    S = SpeciesInteractionNetworks.richness(N)
    l_s = L / S

    A = Matrix(N.edges.edges)

    tls = trophic_level(Bool.(A))

    top = sum(vec(sum(A, dims = 1) .== 0))

    #chain = chain_metrics(A; max_depth = 6)

    # for centrality - freeman centralisation
    c = collect(values(centrality(N)))
    Cmax = maximum(c)

    D = Dict{Symbol,Any}(
        :richness => S,
        :connectance => SpeciesInteractionNetworks.connectance(N),
        :complexity => complexity(N),
        :max_trophic_level => findmax(tls)[1],
        :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
        :generality => std(gen / l_s),
        :vulnerability => std(vul / l_s),
        :top => top / S,
        #:ChLen => chain.ChLen,
        :centrality => sum(Cmax .- c),
        :clustering => clustering(A),
        :trophicCoherence => trophic_coherence(N)
    )

    return D
end

"""
remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    Identifies and sets cannibalistic link to zero
"""
function remove_cannibals(N::SpeciesInteractionNetwork{<:Partiteness,<:Binary})

    # get adj matrix
    S = SpeciesInteractionNetworks.richness(N)
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

function compute_reachable_to_top(N, top_set)

    stack = copy(top)

    reachable[top] .= true

    while !isempty(stack)

        node = pop!(stack)

        for predator in predator_list[node]

            if !reachable[predator]
                reachable[predator] = true
                push!(stack, predator)
            end

        end

    end

    return reachable

end


"""
    chain_metrics(A; max_depth=6)

Calculate food-chain statistics by enumerating all simple food chains from
basal species to top predators.

A food chain is defined as a simple path (no repeated species) beginning
at a basal species and terminating at a top predator.

# Arguments

- `A::AbstractMatrix{Bool}`:
  Adjacency matrix where `A[i,j] == true` indicates predator `i`
  consumes prey `j`.

# Keyword Arguments

- `max_depth=6`:
  Maximum permitted chain length. Recursion terminates once this depth
  is exceeded.

# Returns

Named tuple containing

- `ChLen` : Mean food-chain length.
- `ChSD`  : Standard deviation of food-chain lengths.
- `ChNum` : Logarithm of the number of distinct food chains.

# Notes

This implementation enumerates every simple basal-to-top food chain.
The algorithm is exact but may become computationally expensive for
highly reticulated food webs with many alternative pathways.
"""
function chain_metrics(
    A::AbstractMatrix{Bool};
    max_depth::Int = 6
)

    S = size(A, 1)

    predator_list = [findall(@view A[:, i]) for i in 1:S]

    basal = findall(vec(sum(A; dims = 2) .== 0))
    top = falses(S)
    top[findall(vec(sum(A; dims = 1) .== 0))] .= true

    if isempty(basal) || !any(top)
        return (
            ChLen = NaN,
            ChSD = NaN,
            ChNum = 0.0
        )
    end

    reachable = compute_reachable_to_top(A)

    visited = falses(S)

    count = Ref(0)
    total = Ref(0.0)
    total2 = Ref(0.0)

    function dfs(node::Int, depth::Int)

        if depth > max_depth
            return
        end

        if visited[node]
            return
        end

        if !reachable[node]
            return
        end

        visited[node] = true

        if top[node]

            count[] += 1
            total[] += depth
            total2[] += depth^2

        end

        for predator in predator_list[node]
            dfs(predator, depth + 1)
        end

        visited[node] = false

    end

    for b in basal
        dfs(b, 0)
    end

    if count[] == 0
        return (
            ChLen = NaN,
            ChSD = NaN,
            ChNum = 0.0
        )
    end

    μ = total[] / count[]

    σ = if count[] > 1
        sqrt((total2[] - total[]^2 / count[]) / (count[] - 1))
    else
        0.0
    end

    return (
        ChLen = μ,
        ChSD = σ,
        ChNum = log(count[])
    )

end

"""
trophic_coherence(N::SpeciesInteractionNetwork)

Returns the trophic incoherence parameter q.
Lower q indicates higher trophic coherence.
"""
function trophic_coherence(N::SpeciesInteractionNetwork)

    A = Matrix(N.edges.edges)
    Spp = species(N)
    tl = trophic_level(Bool.(A); species = Spp)

    spp = species(N)
    s = [tl[k] for k in spp]

    trophic_dist = Float64[]

    for i in eachindex(spp)
        for j in eachindex(spp)

            if A[i, j] == true
                push!(trophic_dist, s[i] - s[j])
            end

        end
    end

    # variance of trophic distances
    q = std(trophic_dist)

    return q
end


"""
clustering(A::Matrix{Bool})

    Returns the mean clustering coefficient
"""
function clustering(A::Matrix{Bool})

    N = size(A, 1)
    
    # Calculate the Undirected Degree (k_i)
    # K_i is the total number of neighbors (in-degree + out-degree).
    A_undir = (A + A') .> 0 # A_undir[i,j] = 1 if there is a link i<->j or i->j or i<-j

    # The undirected degree k_i for species i is the sum of the i-th row (or column) of A_undir.
    k_undir = sum(A_undir, dims=2)[:]
    
    # Calculate the Number of Triangles (T_i)
    # In an undirected graph, the number of triangles T_i involving node i is half the (i, i) entry of A_undir^3.
    # We can calculate the total number of undirected links between neighbors of i directly.
    # The element (A_undir^2)_{ij} is the number of 2-paths between i and j.
    # The number of triangles T_i is the sum of links between the neighbors of i.
    
    # Let D be the number of cycles of length 3 (triangles)
    D = diag(A_undir^3) ./ 2

    # Calculate the Local Clustering Coefficient (C_i)
    C_values = Float64[] # Store local clustering coefficients

    for i in 1:N
        k_i = k_undir[i]
        
        # Denominator: Number of possible 2-paths (connections between neighbors)
        # This is the number of pairs of neighbors: k_i * (k_i - 1) / 2
        denominator = k_i * (k_i - 1) / 2
        
        if denominator == 0
            # Species with degree 0 or 1 cannot be part of a triangle.
            push!(C_values, 0.0) 
            continue
        end

        T_i = D[i] # Number of completed triangles involving node i
        
        # Local Clustering Coefficient C_i
        C_i = T_i / denominator
        push!(C_values, C_i)
    end
    
    # Calculate the Mean Clustering Coefficient
    mean_C = mean(C_values)
    
    return mean_C
end

"""
trophic_coherence(N::SpeciesInteractionNetwork)

Returns the trophic incoherence parameter q.
Lower q indicates higher trophic coherence.
"""
function trophic_coherence(N::SpeciesInteractionNetwork)

    A = Matrix(N.edges.edges)
    Spp = species(N)
    tl = trophic_level(Bool.(A); species = Spp)

    spp = species(N)
    s = [tl[k] for k in spp]

    trophic_dist = Float64[]

    for i in eachindex(spp)
        for j in eachindex(spp)

            if A[i, j] == true
                push!(trophic_dist, s[i] - s[j])
            end

        end
    end

    # variance of trophic distances
    q = std(trophic_dist)

    return q
end
