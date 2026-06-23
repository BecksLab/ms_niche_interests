#=
-------------------------------------------------
02_InitialRealism.jl

Assess ecological plausibility of generated food webs.
This script calculates a set of ecological realism
metrics for both the initial and post-dynamics
network states.

Metrics include:
- Species richness (S)
- Connectance
- Proportion of basal species
- Number of isolated species
- Number of illogical consumers
  (species lacking prey but not classified as producers)
- Number of closed feeding loops
- Maximum trophic level
- Mean food-chain length

Outputs:
- realism_initial.csv
- realism_end.csv
-------------------------------------------------
=#

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using LinearAlgebra
using Plots
using SpeciesInteractionNetworks
using Graphs

# --- 2. Load All Code ---
include(joinpath("lib", "internals.jl"))

# --- 3. Custom Metrics ---
function calculate_custom_metrics(
    A::Matrix{Int},
    met_classes::Union{Vector{Symbol}, Missing}
)

    S = size(A, 1)
    L = sum(A)

    conn = S > 0 ? L / (S^2) : 0.0

    in_deg = vec(sum(A, dims = 1))   # prey items per consumer
    out_deg = vec(sum(A, dims = 2))  # predators per species

    isolated = sum((in_deg .== 0) .& (out_deg .== 0))

    prop_basal = S > 0 ? sum(in_deg .== 0) / S : 0.0

    illogical = missing
    if !ismissing(met_classes) && S == length(met_classes)
        illogical = sum(
            (in_deg .== 0) .&
            (met_classes .!= :producer)
        )
    end

    g = SimpleDiGraph(A)
    closed_loops = length(simplecycles(g))

    return (
        S,
        conn,
        prop_basal,
        isolated,
        illogical,
        closed_loops
    )
end

# --- 4. Generic Processor ---
function calculate_realism_metrics(
    input_file::AbstractString,
    output_file::AbstractString
)

    networks = load_object(input_file)

    realism = DataFrame(
        fw_ID = String[],
        model = String[],
        S = Union{Int, Missing}[],
        connectance = Union{Float64, Missing}[],
        prop_basal = Union{Float64, Missing}[],
        mx_tl = Union{Float64, Missing}[],
        mean_FCL = Union{Float64, Missing}[],
        closed_loops = Union{Int, Missing}[],
        isolated = Union{Int, Missing}[],
        illogical = Union{Int, Missing}[]
    )

    for i in 1:nrow(networks)

        fw_id = networks.fw_ID[i]
        model = networks.Model[i]
        A = networks.AdjacencyMatrix[i]

        # Handling of missing networks
        if ismissing(A)

            push!(realism, (
                fw_ID = fw_id,
                model = model,
                S = missing,
                connectance = missing,
                prop_basal = missing,
                mx_tl = missing,
                mean_FCL = missing,
                closed_loops = missing,
                isolated = missing,
                illogical = missing
            ))

            continue
        end

        # build valid network
        net = build_network(A)

        met_classes = networks.MetabolicClasses[i]

        S, conn, p_basal, iso, illog, loops =
            calculate_custom_metrics(A, met_classes)

        # high-level metric summary
        d_tl = missing
        d_fcl = missing

        try
            d = network_summary(net)
            d_tl = d[:max_trophic_level]
            d_fcl = d[:ChLen]
        catch
            @warn "network_summary failed for $(fw_id)"
        end

        # store
        push!(realism, (
            fw_ID = fw_id,
            model = model,
            S = S,
            connectance = conn,
            prop_basal = p_basal,
            mx_tl = d_tl,
            mean_FCL = d_fcl,
            closed_loops = loops,
            isolated = iso,
            illogical = illog
        ))
    end

    CSV.write(output_file, realism)

    return realism
end

# --- 5. Run Analyses ---
realism_initial = calculate_realism_metrics(
    joinpath("networks", "networks.jld2"),
    joinpath("outputs", "realism_initial.csv")
)

realism_end = calculate_realism_metrics(
    joinpath("networks", "networks_END.jld2"),
    joinpath("outputs", "realism_END.csv")
)