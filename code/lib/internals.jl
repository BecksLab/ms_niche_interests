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

    tl = trophic_level(N)

    A = Matrix(N.edges.edges)

    top = sum(vec(sum(A, dims = 1) .== 0))

    chain = chain_metrics(N; max_depth = 6)

    # for centrality - freeman centralisation
    c = collect(values(centrality(N)))
    Cmax = maximum(c)

    D = Dict{Symbol,Any}(
        :richness => S,
        :connectance => SpeciesInteractionNetworks.connectance(N),
        :complexity => complexity(N),
        :trophic_level => mean(collect(values(tl))),
        :distance => distancetobase(N, collect(keys(_gen))[ind_maxgen]),
        :generality => std(gen / l_s),
        :vulnerability => std(vul / l_s),
        :top => top / S,
        :ChLen => chain.ChLen,
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

"""
_trophic_level(N::SpeciesInteractionNetwork)

    Calculates the trophic level of all species in a network using the average 
    shortest path from the prey of species 𝑖 to a basal species

    Williams, Richard J., and Neo D. Martinez. 2004. “Limits to Trophic Levels 
    and Omnivory in Complex Food Webs: Theory and Data.” The American Naturalist 
    163 (3): 458–68. https://doi.org/10.1086/381964.
"""
function trophic_level(N::SpeciesInteractionNetwork)

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

function compute_reachable_to_top(N, top_set)

    reachable = Set(top_set)
    changed = true

    while changed
        changed = false

        for s in species(N)
            if any(n in reachable for n in predecessors(N, s))
                if !(s in reachable)
                    push!(reachable, s)
                    changed = true
                end
            end
        end
    end

    return reachable
end

function chain_metrics(N; max_depth=6)

    # --- STEP 1: Identify basal and top ---
    gen = SpeciesInteractionNetworks.generality(N)
    basal = collect(keys(filter(((k, v),) -> v == 0, gen)))

    vul = SpeciesInteractionNetworks.vulnerability(N)
    top_set = Set(keys(filter(((k, v),) -> v == 0, vul)))

    # If no structure exists → return early
    if isempty(basal) || isempty(top_set)
        return (ChLen = NaN, ChSD = NaN, ChNum = 0.0)
    end

    # --- STEP 2: Reachability pruning ---
    reachable = compute_reachable_to_top(N, top_set)

    # --- STEP 3: Memoized DFS ---
    memo = Dict{Any, Vector{Int}}()

    function dfs(node, visited, depth)

        if depth > max_depth
            return Int[]
        end

        if node in visited
            return Int[]
        end

        # prune unreachable nodes
        if node ∉ reachable
            return Int[]
        end

        if haskey(memo, node)
            return memo[node]
        end

        push!(visited, node)

        lengths = Int[]

        # If this is a top node → chain ends here
        if node in top_set
            push!(lengths, 0)
        end

        # Traverse UP the food web
        for nxt in predecessors(N, node)
            sub_lengths = dfs(nxt, visited, depth + 1)
            for l in sub_lengths
                push!(lengths, l + 1)
            end
        end

        delete!(visited, node)

        memo[node] = lengths
        return lengths
    end

    # --- STEP 4: Collect all chain lengths ---
    all_lengths = Int[]

    for b in basal
        append!(all_lengths, dfs(b, Set(), 0))
    end

    # --- STEP 5: Return summary stats ---
    if isempty(all_lengths)
        return (ChLen = NaN, ChSD = NaN, ChNum = 0.0)
    end

    return (
        ChLen = mean(all_lengths),
        ChSD  = std(all_lengths),
        ChNum = log(length(all_lengths))
    )
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
    tl = trophic_level(N)

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