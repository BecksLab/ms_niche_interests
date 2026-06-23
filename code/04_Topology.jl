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
using Plots
using SpeciesInteractionNetworks

# --- 2. Load All Code ---
include("lib/internals.jl");

# --- 3. functions ---

function calculate_topology(networks::DataFrame)

    topology = DataFrame(
        fw_ID = String[],
        model = String[],
        richness = Union{Int, Missing}[],
        connectance = Union{Float64, Missing}[],
        complexity = Union{Float64, Missing}[],
        max_trophic_level = Union{Float64, Missing}[],
        generality = Union{Float64, Missing}[],
        vulnerability = Union{Float64, Missing}[],
        top = Union{Float64, Missing}[],
        ChLen = Union{Float64, Missing}[],
        distance = Union{Float64, Missing}[],
        centrality = Union{Float64, Missing}[],
        clustering = Union{Float64, Missing}[],
        trophicCoherence = Union{Float64, Missing}[]
    )

    for i in 1:nrow(networks)

        A = networks.AdjacencyMatrix[i]
        model = networks.Model[i]
        id = networks.fw_ID[i]

        # deal with missing networks
        if ismissing(A)

            push!(topology, (
                fw_ID = id,
                model = model,
                richness = missing,
                connectance = missing,
                complexity = missing,
                max_trophic_level = missing,
                generality = missing,
                vulnerability = missing,
                top = missing,
                ChLen = missing,
                distance = missing,
                centrality = missing,
                clustering = missing,
                trophicCoherence = missing
            ))

            continue
        end

        # build network
        net = try
            build_network(A)
        catch
            @warn "build_network failed for model=$model"
            push!(topology, (
                fw_ID = id,
                model = model,
                richness = missing,
                connectance = missing,
                complexity = missing,
                max_trophic_level = missing,
                generality = missing,
                vulnerability = missing,
                top = missing,
                ChLen = missing,
                distance = missing,
                centrality = missing,
                clustering = missing,
                trophicCoherence = missing
            ))
            continue
        end

        # summary
        try
            d = network_summary(net)

            d[:model] = model
            push!(topology, d)

        catch
            @warn "network_summary failed for model=$model"

            push!(topology, (
                fw_ID = id,
                model = model,
                richness = missing,
                connectance = missing,
                complexity = missing,
                max_trophic_level = missing,
                generality = missing,
                vulnerability = missing,
                top = missing,
                ChLen = missing,
                distance = missing,
                centrality = missing,
                clustering = missing,
                trophicCoherence = missing
            ))
        end
    end

    return topology
end

networks = load_object("networks/networks.jld2")
networks_end = load_object("networks/networks_END.jld2")

top_initial = calculate_topology(networks)
top_end = calculate_topology(networks_end)

CSV.write("outputs/topology_initial.csv", top_initial)
CSV.write("outputs/topology_END.csv", top_end)