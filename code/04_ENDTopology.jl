#=
-------------------------------------------------
04_ENDTopology.jl
Calculating structure/topology for the final adj network at equilibrium.
-------------------------------------------------

=#
# Load the Pkg package manager.
using Pkg
# Activate the current project environment to ensure consistent package versions.
Pkg.activate(".")

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using LinearAlgebra
using Plots
using SpeciesInteractionNetworks

# --- 2. Load All Code ---
include("lib/internals.jl");

# --- 3. Import networks .jld2 object ---
path = joinpath(@__DIR__, "data", "outputs", "END_final_adj_03_05_2026.jld2")
networks = load_object(path)

# --- 4. Convert networks to `SpeciesInteractionNetworks` networks ---

# we will also append this to the network object as a col.
networks.InteractionNetwork = [
    ismissing(A) ? missing : build_network(Int.(A))
    for A in networks.alive_connected_A
]

# --- 5. Get topology ---

# create df to store outputs

topology = DataFrame(
    fw_ID = String[],
    model = String[],
    richness = Int64[],
    connectance = Float64[],
    complexity = Float64[],
    trophic_level = Float64[],
    generality = Float64[],
    vulnerability = Float64[],
    top = Float64[],
    ChLen = Float64[],
    distance = Float64[],
    centrality = Float64[],
    clustering = Float64[],
    trophicCoherence = Float64[],
);

# note for simplicity we use a wrapper function

for i in 1:nrow(networks)

    N = networks.InteractionNetwork[i]

    if ismissing(N)
        continue
    end

    try
        d = network_summary(N)

        d[:fw_ID] = networks.fw_ID[i]
        d[:model] = networks.model[i]

        push!(topology, d)

    catch e
        @warn "network_summary failed; skipping row" i fw_ID=networks.fw_ID[i] model=networks.model[i] exception=e
    end

end

# write summaries as .csv

outdir = joinpath(@__DIR__, "data", "outputs")
outpath = joinpath(outdir, "END_topology_03_05_2026.csv")
CSV.write(outpath, topology)
