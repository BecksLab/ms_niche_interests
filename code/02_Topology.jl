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
using DataFrames
using JLD2            # For saving the full data with matrices
using SpeciesInteractionNetworks   # For Uniform, LogUniform distributions

# --- 2. Load All Code ---

# --- 3. Import networks .jld2 object ---

networks = load_object("data/outputs/network_test_verified_seed_42_29-10-2025.jld2")

# --- 4. Convert networks to `SpeciesInteractionNetworks` networks ---

S = networks.S[1]

spp_list = Symbol.([0:1:10])

networks.AdjacencyMatrix[1]