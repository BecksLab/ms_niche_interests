#=
-------------------------------------------------
01_BuildNetworks.jl
Main execution script for generating food webs with the 6 generative models + MaxEnt.
-------------------------------------------------

This script is a "clean" executor. It defines NO functions.
It loads all functions from included files.

It defines all experimental parameters and runs one of two
workflows based on the `verify_web` toggle:

A) `verify_web = false`:
   Runs a `for` loop for `N_REPLICATES` fixed attempts.
   Saves all generated networks, regardless of their properties.
   `run_ID` and `fw_ID` will be identical (e.g., ltm_1, ltm_2...).

B) `verify_web = true`:
   Runs `while` loops that "retry" until `N_REPLICATES`
   valid networks are generated for each model.
   Validity is determined by emergent properties (%basal and connectance) being within the defined filter ranges.
   This range is designed to mimic ecologically "realistic" food webs.
   `run_ID` tracks the attempt number for reproducibility. To be used with the seed.
   `fw_ID` tracks the success number (e.g., ltm_1, ltm_2...).
   A `MAX_ATTEMPTS_MULTIPLIER` is used to prevent infinite loops.
=#
# Load the Pkg package manager.
using Pkg
# Activate the current project environment to ensure consistent package versions.
Pkg.activate(".")

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
using ProgressMeter   # For tracking MaxEnt optimisation
using FoodWebTools    # For network generating models

# --- 2. Load All Code ---

# Get the directory of this script to build relative paths.
SCRIPT_DIR = @__DIR__
# SCRIPT_DIR = "/Users/lauralandonblake/Desktop/Network_buddies/ms_niche_interests/code" # Uncomment for testing if @__DIR__ fails

# Load helper functions
include(joinpath("lib", "verifying_networks.jl"))

# Load MaxEnt model logic
include(joinpath("lib", "maxentmodel.jl"))

# THIS SHOULD NO LONGER BE NEEDED
# Define path to model scripts
#MODEL_SRC_PATH = "/Users/lauralandonblake/Desktop/Network_buddies/FoodWebTools.jl/src/generative_models"

# Load generative model functions
# include(joinpath(MODEL_SRC_PATH, "adbm.jl"))
# include(joinpath(MODEL_SRC_PATH, "cascade_model.jl")) 
# include(joinpath(MODEL_SRC_PATH, "lmatrix.jl"))
# include(joinpath(MODEL_SRC_PATH, "ltm.jl"))
# include(joinpath(MODEL_SRC_PATH, "niche_model.jl"))   
# include(joinpath(MODEL_SRC_PATH, "random_model.jl"))  

# --- 3. Define Global Parameters ---

# Set the master seed for the entire experiment
MASTER_SEED = 42
Random.seed!(MASTER_SEED)

S = 15               # Species richness
N_REPLICATES = 100   # Target number of valid replicates per model

# Output path
OUTPUT_DIR = "data/outputs"

# --- Verification Toggle ---
verify_web = true

# --- Max Attempts Multiplier ---
MAX_ATTEMPTS_MULTIPLIER = 100

# --- Filter Ranges ---
BASAL_RANGE = (0.1, 0.3)        
CONNECTANCE_RANGE = (0.05, 0.3) 

# --- Body Mass Ranges (Log10) ---
BODYMASS_RANGES = (
    producer = (min = -6.0, max = 5.5),
    invertebrate = (min = -5.0, max = 5.7),
)

# --- 4. Initialise Results Storage ---

@info "Initialising master results DataFrame..."

master_df = DataFrame(
    run_ID = Int[],             
    fw_ID = String[],           
    Model = String[],
    S = Int[],
    BasalInput = Union{Float64,Missing}[],     
    ConnectanceInput = Union{Float64,Missing}[], 
    EmergentBasal = Float64[],
    EmergentConnectance = Float64[],
    AdjacencyMatrix = Matrix{Int}[],
    BodyMasses = Union{Vector{Float64},Missing}[],
    MetabolicClasses = Union{Vector{Symbol},Missing}[],
)

# Helper for MaxEnt internal metrics
function calculate_maxent_metrics(adj)
    S_local = size(adj, 1)
    L = sum(adj)
    conn = L / (S_local^2)
    basal_count = sum(sum(adj, dims=1) .== 0)
    p_basal = basal_count / S_local
    return (p_basal = p_basal, conn = conn)
end

# --- 5. Main Execution Block ---

if verify_web
    # --- 5.A. WORKFLOW B: RETRY UNTIL VALID ---
    @info "Starting Workflow B: Retrying until $N_REPLICATES valid networks per model."
    
    MAX_TOTAL_ATTEMPTS_BM = N_REPLICATES * MAX_ATTEMPTS_MULTIPLIER
    MAX_TOTAL_ATTEMPTS_C = N_REPLICATES * MAX_ATTEMPTS_MULTIPLIER

    # --- Loop 1: Body Mass Models ---
    @info "Running Loop 1 (Body Mass Models)..."
    counters_bm = Dict("LTM" => 0, "ATN" => 0, "ADBM" => 0)
    run_ID_bm = 0 
    spp_list_int = 1:S
    spp_list_any = convert(Vector{Any}, spp_list_int) 

    while any(values(counters_bm) .< N_REPLICATES) && run_ID_bm < MAX_TOTAL_ATTEMPTS_BM
        run_ID_bm += 1 
        target_basal_fraction = rand(Uniform(BASAL_RANGE...))
        inputs = generate_bodymass_inputs(S, target_basal_fraction, BODYMASS_RANGES)
        biomass_adbm = inputs.bodymasses .^ (-0.75) 

        result_ltm = run_and_filter_ltm(spp_list_int, inputs.bodymasses, inputs.metabolic_classes, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        result_atn = run_and_filter_atn(spp_list_any, inputs.bodymasses, inputs.is_producer, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        result_adbm = run_and_filter_adbm(spp_list_any, inputs.bodymasses, inputs.is_producer, biomass_adbm, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)

        all_succeeded_bm = (result_ltm !== nothing) && (result_atn !== nothing) && (result_adbm !== nothing)

        if all_succeeded_bm
            if counters_bm["LTM"] < N_REPLICATES
                counters_bm["LTM"] += 1 
                fw_id = "ltm_$(counters_bm["LTM"])" 
                push!(master_df, (run_ID_bm, fw_id, "LTM", S, target_basal_fraction, missing, result_ltm.percent_basal, result_ltm.connectance, result_ltm.adj, inputs.bodymasses, inputs.metabolic_classes))
            end
            if counters_bm["ATN"] < N_REPLICATES
                counters_bm["ATN"] += 1
                fw_id = "atn_$(counters_bm["ATN"])"
                push!(master_df, (run_ID_bm, fw_id, "ATN", S, target_basal_fraction, missing, result_atn.percent_basal, result_atn.connectance, result_atn.adj, inputs.bodymasses, inputs.metabolic_classes))
            end
            if counters_bm["ADBM"] < N_REPLICATES
                counters_bm["ADBM"] += 1
                fw_id = "adbm_$(counters_bm["ADBM"])"
                push!(master_df, (run_ID_bm, fw_id, "ADBM", S, target_basal_fraction, missing, result_adbm.percent_basal, result_adbm.connectance, result_adbm.adj, inputs.bodymasses, inputs.metabolic_classes))
            end
        end
        if run_ID_bm % 500 == 0
            @info "Loop 1 RunID: $run_ID_bm | LTM: $(counters_bm["LTM"]) | ATN: $(counters_bm["ATN"]) | ADBM: $(counters_bm["ADBM"])"
        end
    end

    # --- Loop 2: Connectance Models ---
    @info "Running Loop 2 (Connectance Models)..."
    counters_c = Dict("Niche" => 0, "Cascade" => 0, "Random" => 0)
    run_ID_c = 0 

    while any(values(counters_c) .< N_REPLICATES) && run_ID_c < MAX_TOTAL_ATTEMPTS_C
        run_ID_c += 1
        C_rand = rand(Uniform(CONNECTANCE_RANGE...))

        result_niche = run_and_filter_niche(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        result_cascade = run_and_filter_cascade(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        result_random = run_and_filter_random(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)

        all_succeeded_c = (result_niche !== nothing) && (result_cascade !== nothing) && (result_random !== nothing)

        if all_succeeded_c
            if counters_c["Niche"] < N_REPLICATES
                counters_c["Niche"] += 1
                fw_id = "niche_$(counters_c["Niche"])"
                push!(master_df, (run_ID_c, fw_id, "Niche", S, missing, C_rand, result_niche.percent_basal, result_niche.connectance, result_niche.adj, missing, missing))
            end
            if counters_c["Cascade"] < N_REPLICATES
                counters_c["Cascade"] += 1
                fw_id = "cascade_$(counters_c["Cascade"])"
                push!(master_df, (run_ID_c, fw_id, "Cascade", S, missing, C_rand, result_cascade.percent_basal, result_cascade.connectance, result_cascade.adj, missing, missing))
            end
            if counters_c["Random"] < N_REPLICATES
                counters_c["Random"] += 1
                fw_id = "random_$(counters_c["Random"])"
                push!(master_df, (run_ID_c, fw_id, "Random", S, missing, C_rand, result_random.percent_basal, result_random.connectance, result_random.adj, missing, missing))
            end
        end
        if run_ID_c % 500 == 0
            @info "Loop 2 RunID: $run_ID_c | Niche: $(counters_c["Niche"]) | Cascade: $(counters_c["Cascade"]) | Random: $(counters_c["Random"])"
        end
    end

    # --- Loop 3: MaxEnt Model ---
    @info "Running Loop 3 (MaxEnt Model)..."
    success_count_me = 0
    run_ID_me = 0
    p_me = Progress(N_REPLICATES, dt=0.5, desc="MaxEnt Gen: ", color=:cyan)

    while success_count_me < N_REPLICATES && run_ID_me < (N_REPLICATES * MAX_ATTEMPTS_MULTIPLIER)
        run_ID_me += 1
        C_target = rand(Uniform(CONNECTANCE_RANGE...))
        L_target = round(Int, C_target * S^2)
        
        jdd = joint_degree_dist_maxent(S, L_target)
        local degrees
        valid_seq = false
        for _ in 1:100 
            degrees = simulate_degrees(jdd)
            if sum(degrees.kin) == sum(degrees.kout) && sum(degrees.kin) > 0
                valid_seq = true
                break
            end
        end
        !valid_seq && continue 

        adj = build_maxent_network(degrees.kin, degrees.kout; iterations=5000, T0=0.1)
        metrics = calculate_maxent_metrics(adj)
        
        if metrics.p_basal >= BASAL_RANGE[1] && metrics.p_basal <= BASAL_RANGE[2] &&
           metrics.conn >= CONNECTANCE_RANGE[1] && metrics.conn <= CONNECTANCE_RANGE[2]
            
            success_count_me += 1
            push!(master_df, (run_ID_me, "maxent_$(success_count_me)", "MaxEnt", S, missing, C_target, metrics.p_basal, metrics.conn, Int.(adj), missing, missing))
            next!(p_me; showvalues = [(:Attempts, run_ID_me), (:Successes, success_count_me)])
        end
    end

else
    # --- 5.B. WORKFLOW A: FIXED ATTEMPTS ---
    @info "Starting Workflow A: Performing $N_REPLICATES fixed attempts per model. No filtering."
    
    spp_list_int = 1:S
    spp_list_any = convert(Vector{Any}, spp_list_int)

    for i = 1:N_REPLICATES
        target_basal_fraction = rand(Uniform(BASAL_RANGE...))
        inputs = generate_bodymass_inputs(S, target_basal_fraction, BODYMASS_RANGES)
        biomass_adbm = inputs.bodymasses .^ (-0.75)
        C_rand = rand(Uniform(CONNECTANCE_RANGE...))
        run_ID = i 

        # BM Models
        result_ltm = run_and_filter_ltm(spp_list_int, inputs.bodymasses, inputs.metabolic_classes, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "ltm_$i", "LTM", S, target_basal_fraction, missing, result_ltm.percent_basal, result_ltm.connectance, result_ltm.adj, inputs.bodymasses, inputs.metabolic_classes))
        
        result_atn = run_and_filter_atn(spp_list_any, inputs.bodymasses, inputs.is_producer, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "atn_$i", "ATN", S, target_basal_fraction, missing, result_atn.percent_basal, result_atn.connectance, result_atn.adj, inputs.bodymasses, inputs.metabolic_classes))

        result_adbm = run_and_filter_adbm(spp_list_any, inputs.bodymasses, inputs.is_producer, biomass_adbm, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "adbm_$i", "ADBM", S, target_basal_fraction, missing, result_adbm.percent_basal, result_adbm.connectance, result_adbm.adj, inputs.bodymasses, inputs.metabolic_classes))

        # Connectance Models
        result_niche = run_and_filter_niche(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "niche_$i", "Niche", S, missing, C_rand, result_niche.percent_basal, result_niche.connectance, result_niche.adj, missing, missing))
        
        result_cascade = run_and_filter_cascade(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "cascade_$i", "Cascade", S, missing, C_rand, result_cascade.percent_basal, result_cascade.connectance, result_cascade.adj, missing, missing))

        result_random = run_and_filter_random(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "random_$i", "Random", S, missing, C_rand, result_random.percent_basal, result_random.connectance, result_random.adj, missing, missing))

        # MaxEnt Models
        L_target = round(Int, C_rand * S^2)
        jdd = joint_degree_dist_maxent(S, L_target)
        local deg_samp
        while true # Ensure we get a valid sequence even in Workflow A
            deg_samp = simulate_degrees(jdd)
            if sum(deg_samp.kin) == sum(deg_samp.kout) && sum(deg_samp.kin) > 0; break; end
        end
        adj_me = build_maxent_network(deg_samp.kin, deg_samp.kout; iterations=3000)
        m_me = calculate_maxent_metrics(adj_me)
        push!(master_df, (run_ID, "maxent_$i", "MaxEnt", S, missing, C_rand, m_me.p_basal, m_me.conn, Int.(adj_me), missing, missing))
    end
end 

# --- 6. Save Final Results ---

@info "All replicates generated. Saving results..."

status_str = verify_web ? "verified" : "unverified"
date_str = Dates.format(now(), "dd-mm-yyyy")
filename_base = "network_test_$(status_str)_seed_$(MASTER_SEED)_$(date_str)"

mkpath(OUTPUT_DIR)
csv_path = joinpath(OUTPUT_DIR, "$(filename_base).csv")
jld2_path = joinpath(OUTPUT_DIR, "$(filename_base).jld2")

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

@info "--- Workflow Complete ---"