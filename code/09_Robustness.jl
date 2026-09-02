#=
-------------------------------------------------
09_Robustness.jl
Quick interlude on R50 values.
-------------------------------------------------

=#

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using Extinctions

# --- 2. Load All Code ---
include("lib/internals.jl");

# load networks
networks = load_object("networks/networks.jld2")

# calculate robustness

robustness_vals = DataFrame(
        Model = String[],
        fw_ID = String[],
        run_ID = Int[],
        robustness = Float64[],
    )

@showprogress "Calculating robustness" for i in 1:nrow(networks)

    N = build_network(networks.AdjacencyMatrix[i])
    N = remove_cannibals(N)

    rob = robustness(
        N;
        threshold = 50,
        mechanism = :cascade,
    )

    push!(robustness_vals, (
        Model = networks.Model[i],
        fw_ID = networks.fw_ID[i],
        run_ID = networks.run_ID[i],
        robustness = rob,
    ))
end

CSV.write("outputs/robustness.csv", robustness_vals)
