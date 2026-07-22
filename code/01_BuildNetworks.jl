#=
-------------------------------------------------
01_BuildNetworks.jl
Main execution script for generating food webs with the 6 generative models + MaxEnt.
-------------------------------------------------

This script generates exactly N_REPLICATES webs per model without applying 
topological filters to the emergent properties. This allows us to observe 
the true structural distribution (including illogical/isolated species) 
produced by each generative algorithm.
=#

# --- 1. Load Dependencies ---
using DataFrames      # For managing results in a table
using CSV             # For saving the summary CSV file
using JLD2            # For saving the full data with matrices
using Distributions   # For Uniform, LogUniform distributions
using Random          # For setting the seed
using Dates           # For timestamping the output file
using Statistics      # Required by helpers
using LinearAlgebra   # Required by helpers
using Graphs          # Required by helpers/models
using Ipopt           # Required by MaxEnt model optimisation
import JuMP: Model    # Required by MaxEnt model optimisation
using StatsBase       # For sampling functions
using ProgressMeter   # For tracking MaxEnt optimisation
using FoodWebTools    # For network generating models

# --- 2. Load All Code ---

# Load helper functions using the absolute path directory
include(joinpath("lib", "verifying_networks.jl"))

# Load MaxEnt model logic using the absolute path directory
include(joinpath("lib", "maxentmodel.jl"))

# --- 3. Define Global Parameters ---

# Set the master seed for the entire experiment
MASTER_SEED = 66
Random.seed!(MASTER_SEED)

SPECIES_RICHNESS = [10, 20, 40]
N_REPLICATES = 100 # number networks per S

# --- Verification Toggle ---
# Set to false to entirely bypass the filtering logic in `lib/verifying_networks.jl`
verify_web = false

# --- Input Parameter Ranges ---
# Note: These are now ONLY used as input parameters to seed the models, 
# not as emergent filters to delete networks.
BASAL_RANGE = (0.1, 0.3)
CONNECTANCE_RANGE = (0.05, 0.3)

# --- Body Mass Ranges (Log10) ---
BODYMASS_RANGES = (
    producer=(min=-6.0, max=5.5),
    invertebrate=(min=-5.0, max=5.7),
)

# --- 4. Initialise Results Storage ---

@info "Initialising master results DataFrame..."

master_df = DataFrame(
    run_ID=Int[],
    fw_ID=String[],
    Model=String[],
    S=Int[],
    BasalInput=Union{Float64,Missing}[],
    ConnectanceInput=Union{Float64,Missing}[],
    EmergentBasal=Float64[],
    EmergentConnectance=Float64[],
    AdjacencyMatrix=Matrix{Int}[],
    BodyMasses=Union{Vector{Float64},Missing}[],
    MetabolicClasses=Union{Vector{Symbol},Missing}[],
)

# Helper for MaxEnt internal metrics
function calculate_maxent_metrics(adj)
    S_local = size(adj, 1)
    L = sum(adj)
    conn = L / (S_local^2)
    basal_count = sum(sum(adj, dims=1) .== 0)
    p_basal = basal_count / S_local
    return (p_basal=p_basal, conn=conn)
end

# --- 5. Main Execution Block ---

@info "Starting Generative Workflow..."

for S in SPECIES_RICHNESS

    @info "Beginning simulations for S = $S"

    spp_list_int = 1:S
    spp_list_any = convert(Vector{Any}, spp_list_int)

    for i = 1:N_REPLICATES

        # 1. Sample inputs for this replicate
        target_basal_fraction = rand(Uniform(BASAL_RANGE...))
        inputs = generate_bodymass_inputs(S, target_basal_fraction, BODYMASS_RANGES)
        biomass_adbm = inputs.bodymasses .^ (-0.75)
        C_rand = rand(Uniform(CONNECTANCE_RANGE...))
        run_ID = i

        # 2. Body Mass Models
        result_ltm = run_and_filter_ltm(spp_list_int,
            inputs.bodymasses,
            inputs.metabolic_classes,
            verify_web,
            BASAL_RANGE,
            CONNECTANCE_RANGE)
        push!(master_df, (run_ID,
            "ltm_$(i)_$(S)",
            "LTM",
            S,
            target_basal_fraction, missing,
            result_ltm.percent_basal,
            result_ltm.connectance,
            result_ltm.adj,
            inputs.bodymasses,
            inputs.metabolic_classes))

        result_atn = run_and_filter_atn(spp_list_any,
            inputs.bodymasses,
            inputs.is_producer,
            verify_web,
            BASAL_RANGE,
            CONNECTANCE_RANGE)
        push!(master_df, (run_ID,
            "atn_$(i)_$(S)",
            "ATN",
            S,
            target_basal_fraction, missing,
            result_atn.percent_basal,
            result_atn.connectance,
            result_atn.adj,
            inputs.bodymasses,
            inputs.metabolic_classes))

        result_adbm = run_and_filter_adbm(spp_list_any,
            inputs.bodymasses,
            inputs.is_producer,
            biomass_adbm,
            verify_web,
            BASAL_RANGE,
            CONNECTANCE_RANGE)
        push!(master_df, (run_ID,
            "adbm_$(i)_$(S)",
            "ADBM",
            S,
            target_basal_fraction,
            missing,
            result_adbm.percent_basal,
            result_adbm.connectance,
            result_adbm.adj,
            inputs.bodymasses,
            inputs.metabolic_classes))

        # 3. Connectance Models
        result_niche = run_and_filter_niche(S,
            C_rand,
            verify_web,
            BASAL_RANGE,
            CONNECTANCE_RANGE)
        push!(master_df, (run_ID,
            "niche_$(i)_$(S)",
            "Niche",
            S,
            missing,
            C_rand,
            result_niche.percent_basal,
            result_niche.connectance,
            result_niche.adj,
            missing,
            missing))

        result_cascade = run_and_filter_cascade(S,
            C_rand,
            verify_web,
            BASAL_RANGE,
            CONNECTANCE_RANGE)
        push!(master_df, (run_ID,
            "cascade_$(i)_$(S)",
            "Cascade",
            S,
            missing,
            C_rand,
            result_cascade.percent_basal,
            result_cascade.connectance,
            result_cascade.adj,
            missing,
            missing))

        result_random = run_and_filter_random(S,
            C_rand,
            verify_web,
            BASAL_RANGE,
            CONNECTANCE_RANGE)
        push!(master_df, (run_ID,
            "random_$(i)_$(S)",
            "Random",
            S,
            missing,
            C_rand,
            result_random.percent_basal,
            result_random.connectance,
            result_random.adj,
            missing,
            missing))

        # 4. MaxEnt Models
        L_target = round(Int, C_rand * S^2)
        jdd = joint_degree_dist_maxent(S, L_target)
        local deg_samp

        # We must maintain this while loop strictly to satisfy the structural rules 
        # of building a matrix (in-degrees == out-degrees), otherwise the builder fails. 
        # This does not filter emergent properties.
        while true
            deg_samp = simulate_degrees(jdd)
            if sum(deg_samp.kin) == sum(deg_samp.kout) && sum(deg_samp.kin) > 0
                ;
                break;
            end
        end
        # The builder can still fail if the degree sequence is not graphical, so we wrap in a try-catch to ensure we get exactly 
        # N_REPLICATES successful builds.
        adj_me = build_maxent_network(deg_samp.kin, deg_samp.kout; iterations=3000)
        # Calculate the emergent metrics for this MaxEnt network to store alongside the summary
        m_me = calculate_maxent_metrics(adj_me)
        push!(master_df, (run_ID, "maxent_$(i)_$(S)",
            "MaxEnt",
            S,
            missing,
            C_rand,
            m_me.p_basal,
            m_me.conn,
            Int.(adj_me),
            missing,
            missing))

        if i % 10 == 0
            @info "Generated $i / $N_REPLICATES replicates across all models."
        end
    end

    @info "Finished S = $S"

end

# --- 6. Save Final Results ---

@info "All replicates generated. Saving un-filtered results..."

csv_path = joinpath("networks", "networks.csv")
jld2_path = joinpath("networks", "networks.jld2")

@info "Saving full data to CSV..."
try
    CSV.write(csv_path, master_df)
    @info "Successfully saved full data to: $csv_path"
catch e
    @warn "Failed to save CSV file!"
    println(e)
end

@info "Saving full data to JLD2..."
try
    JLD2.save_object(jld2_path, master_df)
    @info "Successfully saved full data to: $jld2_path"
catch e
    @warn "Failed to save JLD2 file!"
    println(e)
end

@info "--- Generation Complete ---"