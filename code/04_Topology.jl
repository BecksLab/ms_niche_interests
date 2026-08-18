#=
-------------------------------------------------
04_Topology.jl
Calculating structure/topology for generated networks.
-------------------------------------------------

=#

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using LinearAlgebra
using SpeciesInteractionNetworks
using FoodWebTools

# --- 2. Load All Code ---
include("lib/internals.jl");

# --- 3. functions ---

function calculate_topology(networks::DataFrame)

    topology = DataFrame(
        fw_ID=String[],
        model=String[],
        richness=Union{Int,Missing}[],
        links=Union{Any,Missing}[],
        connectance=Union{Any,Missing}[],
        distance=Union{Any,Missing}[],
        basal=Union{Any,Missing}[],
        top=Union{Any,Missing}[],
        intermediate=Union{Any,Missing}[],
        diameter=Union{Any,Missing}[],
        herbivory=Union{Any,Missing}[],
        omnivory=Union{Any,Missing}[],
        predpreyRatio=Union{Any,Missing}[],
        l_S=Union{Any,Missing}[],
        GenSD=Union{Any,Missing}[],
        VulSD=Union{Any,Missing}[],
        TL=Union{Any,Missing}[],
        ChLen=Union{Any,Missing}[],
        ChSD=Union{Any,Missing}[],
        ChNum=Union{Any,Missing}[],
        path=Union{Any,Missing}[],
        LinkSD=Union{Any,Missing}[],
        S1=Union{Any,Missing}[],
        S2=Union{Any,Missing}[],
        S4=Union{Any,Missing}[],
        S5=Union{Any,Missing}[],
        centrality=Union{Any,Missing}[],
        loops=Union{Any,Missing}[],
        intervals=Union{Any,Missing}[],
        MaxSim=Union{Any,Missing}[],
        Clust=Union{Any,Missing}[],
        trophicCoherence=Union{Any,Missing}[],
        trophicVar=Union{Any,Missing}[],
    )

    for i in 1:nrow(networks)

        A = networks.AdjacencyMatrix[i]
        model = networks.Model[i]
        id = networks.fw_ID[i]

        if ismissing(A)

            push!(topology, (
                    fw_ID=id,
                    model=model,
                    richness=0,
                    links=missing,
                    connectance=missing,
                    distance=missing,
                    basal=missing,
                    top=missing,
                    intermediate=missing,
                    diameter=missing,
                    herbivory=missing,
                    omnivory=missing,
                    predpreyRatio=missing,
                    l_S=missing,
                    GenSD=missing,
                    VulSD=missing,
                    TL=missing,
                    ChLen=missing,
                    ChSD=missing,
                    ChNum=missing,
                    path=missing,
                    LinkSD=missing,
                    S1=missing,
                    S2=missing,
                    S4=missing,
                    S5=missing,
                    centrality=missing,
                    loops=missing,
                    intervals=missing,
                    MaxSim=missing,
                    Clust=missing,
                    trophicCoherence=missing,
                    trophicVar=missing
            ))

        else
            net = build_network(A)

            # guard against empty / degenerate networks
            if isempty(A) || count(!iszero, A) == 0
                push!(topology, (
                    fw_ID=id,
                    model=model,
                    richness=0,
                    links=missing,
                    connectance=missing,
                    distance=missing,
                    basal=missing,
                    top=missing,
                    intermediate=missing,
                    diameter=missing,
                    herbivory=missing,
                    omnivory=missing,
                    predpreyRatio=missing,
                    l_S=missing,
                    GenSD=missing,
                    VulSD=missing,
                    TL=missing,
                    ChLen=missing,
                    ChSD=missing,
                    ChNum=missing,
                    path=missing,
                    LinkSD=missing,
                    S1=missing,
                    S2=missing,
                    S4=missing,
                    S5=missing,
                    centrality=missing,
                    loops=missing,
                    intervals=missing,
                    MaxSim=missing,
                    Clust=missing,
                    trophicCoherence=missing,
                    trophicVar=missing
                ))
                continue
            end

            d = network_summary(net)

            d[:model] = model
            d[:fw_ID] = id
            push!(topology, d)
        end
    end

    return topology
end

networks = load_object("networks/networks.jld2")
networks.AdjacencyMatrix = map(networks.AdjacencyMatrix) do x
    if ismissing(x)
        missing
    elseif isa(x, AbstractMatrix) && isempty(x)
        missing
    else
        Int.(x)
    end
end
networks_end = load_object("networks/networks_END.jld2")
networks_end.AdjacencyMatrix = map(networks_end.AdjacencyMatrix) do x
    if ismissing(x)
        missing
    elseif isa(x, AbstractMatrix) && isempty(x)
        missing
    else
        Int.(x)
    end
end

top_initial = calculate_topology(networks)
top_end = calculate_topology(networks_end)

CSV.write("outputs/topology_initial.csv", top_initial)
CSV.write("outputs/topology_END.csv", top_end)