#=
-------------------------------------------------
04_ENDTopology.jl
Calculating structure/topology for post-simulation networks.
Input: post_networks.jld2
Output: post_topology.csv
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
path = joinpath(@__DIR__, "data", "outputs", "post_networks.jld2")
post_networks = load_object(path)

# --- 4. Convert networks to `SpeciesInteractionNetworks` networks ---

# we will also append this to the network object as a col.
post_networks.InteractionNetwork = [
    ismissing(A) ? missing : build_network(Int.(A))
    for A in post_networks.post_adj
]

# --- 5. Get post_topology ---

# create df to store outputs

post_topology = DataFrame(
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

for i in 1:nrow(post_networks)

    N = post_networks.InteractionNetwork[i]

    if ismissing(N)
        continue
    end

    try
        d = network_summary(N)

        d[:fw_ID] = post_networks.fw_ID[i]
        d[:model] = post_networks.model[i]

        push!(post_topology, d)

    catch e
        @warn "network_summary failed; skipping row" i fw_ID=post_networks.fw_ID[i] model=post_networks.model[i] exception=e
    end

end

# write summaries as .csv

outdir = joinpath(@__DIR__, "data", "outputs")
outpath = joinpath(outdir, "post_topology.csv")
CSV.write(outpath, post_topology)
