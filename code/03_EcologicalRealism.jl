#=
-------------------------------------------------
02_InitialRealism.jl

Assess ecological plausibility of generated food webs.
This script calculates a set of ecological realism
metrics for both the initial and post-dynamics
network states.
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
using FoodWebTools

# --- 2. Load All Code ---
include(joinpath("lib", "internals.jl"))

# --- 3. Custom Metrics ---
function calculate_custom_metrics(
    A::AbstractMatrix{<:Integer},
    met_classes::Union{Vector,Missing};
    max_species_for_loops::Int=15 # Safety check threshold
)

    S = size(A, 1)
    L = sum(A)

    conn = S > 0 ? L / (S^2) : 0.0

    out_deg = vec(sum(A, dims=1))   # pred items per consumer
    in_deg = vec(sum(A, dims=2))  # prey per species

    isolated = sum((in_deg .== 0) .& (out_deg .== 0))

    illogical = missing
    if !ismissing(met_classes) && S == length(met_classes)
        illogical = sum(
            (in_deg .== 0) .&
            (met_classes .!= :producer)
        )
    end

    # 💀: skip cycle counting if network is too large/dense
    # Exponential algorithm will crash the session otherwise.
    closed_loops = missing
    if S <= max_species_for_loops
        try
            g = SimpleDiGraph(A)
            closed_loops = length(simplecycles(g))
        catch e
            @warn "Cycle detection failed or ran out of memory."
            closed_loops = missing
        end
    else
        @warn "Network size ($S species) exceeds safety limit. Skipping closed loop count."
    end

    return (
        S,
        conn,
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
    total_networks = nrow(networks)

    realism = DataFrame(
        fw_ID=String[],
        model=String[],
        S=Union{Int,Missing}[],
        connectance=Union{Float64,Missing}[],
        mx_tl=Union{Float64,Missing}[],
        intervality=Union{Int,Missing}[],
        isolated=Union{Int,Missing}[],
        illogical=Union{Int,Missing}[],
        loops=Union{Int,Missing}[]
    )

    for i in 1:total_networks

        fw_id = networks.fw_ID[i]
        model = networks.Model[i]
        A = networks.AdjacencyMatrix[i]

        # Handling of missing networks
        if ismissing(A)
            push!(realism, (
                fw_ID=fw_id,
                model=model,
                S=missing,
                connectance=missing,
                mx_tl=missing,
                intervality=missing,
                isolated=missing,
                illogical=missing,
                loops=missing
            ))
            continue
        end

        # Build valid network
        A = Int.(A)
        met_classes = networks.MetabolicClasses[i]

        ## Skip failed networks stored as sparse matrices
        if A isa SparseArrays.SparseMatrixCSC

            @warn "Skipping malformed network $fw_id"

            push!(realism, (
                fw_ID=fw_id,
                model=model,
                S=missing,
                connectance=missing,
                mx_tl=missing,
                intervality=missing,
                isolated=missing,
                illogical=missing,
                loops=missing
            ))

            continue
        end

        # S, conn, iso, illog, loops
        S, conn, iso, illog, loop =
            calculate_custom_metrics(A, met_classes; max_species_for_loops=100)

        # high-level metric summary
        d_tl = missing
        interval = missing

        try
            d_tl = findmax(trophic_level(Bool.(A)))[1]
            interval = intervality(Bool.(A))
        catch
            @warn "network_summary failed for $(fw_id)"
        end

        # store
        push!(realism, (
            fw_ID=fw_id,
            model=model,
            S=S,
            connectance=conn,
            mx_tl=d_tl,
            intervality=interval,
            isolated=iso,
            illogical=illog,
            loops=loop
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