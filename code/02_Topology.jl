#=
-------------------------------------------------
02_Topology.jl
Calculating structure/topology for generated networks.
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
using SpeciesInteractionNetworks

# --- 2. Load All Code ---
include("lib/internals.jl");

# --- 3. Import networks .jld2 object ---

networks = load_object("data/outputs/network_test_verified_seed_42_29-10-2025.jld2")

# --- 4. Convert networks to `SpeciesInteractionNetworks` networks ---

# we will also append this to the network object as a col.
networks.InteractionNetwork = build_network.(networks.AdjacencyMatrix)

# --- 5. Get topology ---

# create df to store outputs

topology = DataFrame(
    model = String[],
    richness = Int64[],
    connectance = Float64[],
    diameter = Int64[],
    complexity = Float64[],
    trophic_level = Float64[],
    distance = Float64[],
    generality = Float64[],
    vulnerability = Float64[],
    redundancy = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
);

# note for simplicity we use a wrapper function

for i in 1:nrow(networks)

    d = network_summary(networks.InteractionNetwork[i])

    d[:model] = networks.Model[i]

    push!(topology, d)

end

# write summaries as .csv
CSV.write("data/outputs/topology.csv", topology)