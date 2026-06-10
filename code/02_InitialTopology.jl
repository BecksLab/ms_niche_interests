#=
-------------------------------------------------
02_InitialTopology.jl
Calculating initial structure/topology for generated networks.
-------------------------------------------------
=#

# Load the Pkg package manager.
using Pkg
Pkg.activate(".")

# Explicitly pull SpeciesInteractionNetworks from its GitHub repository
Pkg.add(url="https://github.com/PoisotLab/SpeciesInteractionNetworks.jl.git")

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using LinearAlgebra
using Plots
using SpeciesInteractionNetworks
using Graphs # Required to calculate closed loops (simplecycles)

# --- 2. Load All Code ---

include(joinpath("lib", "internals.jl"));

# --- 3. Import networks .jld2 object ---
# Read the unfiltered output from script 1
date_str = "05-06-2026"
networks = load_object(joinpath(BASE_DIR, "data", "outputs", "network_test_unfiltered_seed_42_$(date_str).jld2"))

# --- 4. Convert networks to `SpeciesInteractionNetworks` ---
networks.InteractionNetwork = build_network.(networks.AdjacencyMatrix)

# --- 5. Custom Topology Function ---
function calculate_initial_custom_metrics(A::Matrix{Int}, met_classes::Union{Vector{Symbol}, Missing})
    S = size(A, 1)
    L = sum(A)
    conn = S > 0 ? L / (S^2) : 0.0
    
    in_deg = vec(sum(A, dims=1))  # Prey items per consumer
    out_deg = vec(sum(A, dims=2)) # Predators per species
    
    # 1. Completely isolated nodes (0 links)
    isolated = sum((in_deg .== 0) .& (out_deg .== 0))
    
    # 2. Proportion basal (no prey)
    prop_basal = S > 0 ? sum(in_deg .== 0) / S : 0.0
    
    # 3. Illogical consumers (no prey, but not classified as producers)
    illogical = missing
    if !ismissing(met_classes) && S == length(met_classes)
        illogical = sum((in_deg .== 0) .& (met_classes .!= :producer))
    end
    
    # 4. Closed loops
    g = SimpleDiGraph(A)
    cycles = simplecycles(g)
    closed_loops = length(cycles)
    
    return (S, conn, prop_basal, isolated, illogical, closed_loops)
end

# --- 6. Get topology ---
topology = DataFrame(
    fw_ID = String[],
    model = String[],
    S_initial = Int[],
    connectance_initial = Float64[],
    prop_basal_initial = Float64[],
    mx_tl_initial = Float64[],
    mean_FCL_initial = Float64[],
    closed_loops_initial = Int[],
    isolated_initial = Int[],
    illogical_initial = Union{Int, Missing}[],
);

for i in 1:nrow(networks)
    
    met_classes = ismissing(networks.MetabolicClasses[i]) ? missing : networks.MetabolicClasses[i]
    S_init, conn_init, p_basal_init, iso_init, illog_init, loops_init = calculate_initial_custom_metrics(networks.AdjacencyMatrix[i], met_classes)

    # Some completely disjointed networks may fail default summary metrics
    d_tl = NaN
    d_fcl = NaN
    try
        d = network_summary(networks.InteractionNetwork[i])
        d_tl = d[:trophic_level]
        d_fcl = d[:ChLen]
    catch
        @warn "Standard network_summary failed for $(networks.fw_ID[i]). Assigning NaN to TL and FCL."
    end
    
    push!(topology, (
        fw_ID = networks.fw_ID[i],
        model = networks.Model[i],
        S_initial = S_init,
        connectance_initial = conn_init,
        prop_basal_initial = p_basal_init,
        mx_tl_initial = d_tl, 
        mean_FCL_initial = d_fcl,
        closed_loops_initial = loops_init,
        isolated_initial = iso_init,
        illogical_initial = illog_init
    ))
end

# write summaries as .csv securely using absolute paths
CSV.write(joinpath("data", "outputs", "topology_initial_$(date_str).csv"), topology)