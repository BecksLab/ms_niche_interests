#=
-------------------------------------------------
01_BuildNetworks.jl
Main execution script for generating food webs with the 6 generative models.
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

# --- 2. Load All Code ---

# Get the directory of this script to build relative paths.
# `@__DIR__` is a special macro that gets the path of the current file.
SCRIPT_DIR = @__DIR__
# SCRIPT_DIR = "/Users/lauralandonblake/Desktop/Network_buddies/ms_niche_interests/code" # Uncomment for testing if @__DIR__ fails

# Load all helper functions, wrappers, and calculators
# This file contains generate_bodymass_inputs, _process_web, run_and_filter_... etc.
include(joinpath(SCRIPT_DIR, "lib", "verifying_networks.jl"))

# Define path to model scripts (assuming a fixed relative structure)
MODEL_SRC_PATH = "/Users/lauralandonblake/Desktop/Network_buddies/FoodWebTools.jl/src/generative_models"

# Load all 6 (pure) generative model functions
# This assumes cascade, niche, and random have been modified to remove internal filtering
include(joinpath(MODEL_SRC_PATH, "adbm.jl"))
include(joinpath(MODEL_SRC_PATH, "cascade_model.jl")) 
include(joinpath(MODEL_SRC_PATH, "lmatrix.jl"))
include(joinpath(MODEL_SRC_PATH, "ltm.jl"))
include(joinpath(MODEL_SRC_PATH, "niche_model.jl"))   
include(joinpath(MODEL_SRC_PATH, "random_model.jl"))  


# --- 3. Define Global Parameters ---

# Set the master seed for the entire experiment for 100% reproducibility.
MASTER_SEED = 42
Random.seed!(MASTER_SEED)

S = 15               # Species richness
N_REPLICATES = 100   # Target number of valid replicates per model

# --- *** Set explicit output path *** ---
OUTPUT_DIR = "/Users/lauralandonblake/Desktop/Network_buddies/ms_niche_interests/code/data/outputs"

# --- Verification Toggle ---
# true = Workflow B: "Retry" loop. Only save networks that pass filters.
# false = Workflow A: "Fixed attempts" loop. Save all generated networks.
verify_web = true

# --- Max Attempts Multiplier ---
# (Only used if verify_web = true)
# The total number of attempts will be N_REPLICATES * MAX_ATTEMPTS_MULTIPLIER
MAX_ATTEMPTS_MULTIPLIER = 100

# --- Filter Ranges ---
# These tuples define the (min, max) allowed values.
BASAL_RANGE = (0.1, 0.3)        # Allowed emergent % basal (min, max)
CONNECTANCE_RANGE = (0.05, 0.3) # Allowed emergent connectance (min, max)

# --- Body Mass Ranges (Log10) ---
# These are the log10 exponents for the body mass ranges.
BODYMASS_RANGES = (
    producer = (min = -6.0, max = 5.5),
    invertebrate = (min = -5.0, max = 5.7),
)

# --- 4. Initialise Results Storage ---

@info "Initialising master results DataFrame..."

# Define the schema for the master DataFrame.
# This structure will hold all results from all models.
master_df = DataFrame(
    run_ID = Int[],             # The attempt number (for reproducibility)
    fw_ID = String[],           # The success ID (e.g., "ltm_1")
    Model = String[],
    S = Int[],
    BasalInput = Union{Float64,Missing}[],     # TargetBasal for BM models
    ConnectanceInput = Union{Float64,Missing}[], # TargetC for C models
    EmergentBasal = Float64[],
    EmergentConnectance = Float64[],
    AdjacencyMatrix = Matrix{Int}[],
    BodyMasses = Union{Vector{Float64},Missing}[],
    MetabolicClasses = Union{Vector{Symbol},Missing}[],
)


# --- 5. Main Execution Block ---

# This block chooses which workflow to run based on the `verify_web` toggle.
if verify_web
    # --- 5.A. WORKFLOW B: RETRY UNTIL VALID (While Loops) ---
    @info "Starting Workflow B: Retrying until $N_REPLICATES valid networks per model."
    
    # Calculate the maximum number of attempts for each loop.
    MAX_TOTAL_ATTEMPTS_BM = N_REPLICATES * MAX_ATTEMPTS_MULTIPLIER
    MAX_TOTAL_ATTEMPTS_C = N_REPLICATES * MAX_ATTEMPTS_MULTIPLIER

    # --- Loop 1: Body Mass Models ---
    @info "Running Loop 1 (Body Mass Models)..."
    # counters_bm tracks the `fw_ID` (successes) for each model.
    counters_bm = Dict("LTM" => 0, "ATN" => 0, "ADBM" => 0)
    # run_ID_bm tracks the `run_ID` (attempts).
    run_ID_bm = 0 
    spp_list_int = 1:S
    spp_list_any = convert(Vector{Any}, spp_list_int) # For adbm/lmatrix

    # Continue looping as long as any model has fewer successes than N_REPLICATES
    # AND we have not exceeded the maximum number of attempts.
    while any(values(counters_bm) .< N_REPLICATES) && run_ID_bm < MAX_TOTAL_ATTEMPTS_BM
        run_ID_bm += 1 # Increment attempt counter
        
        # 1. Generate new shared inputs for this attempt.
        target_basal_fraction = rand(Uniform(BASAL_RANGE...))
        inputs = generate_bodymass_inputs(S, target_basal_fraction, BODYMASS_RANGES)
        biomass_adbm = inputs.bodymasses .^ (-0.75) # ADBM-specific input

        # --- *** NEW LOGIC: "ALL OR NOTHING" *** ---
        # 2. Run all 3 models in the group FIRST.
        #    Do not check counters here.
        result_ltm = run_and_filter_ltm(spp_list_int, inputs.bodymasses, inputs.metabolic_classes, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        result_atn = run_and_filter_atn(spp_list_any, inputs.bodymasses, inputs.is_producer, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        result_adbm = run_and_filter_adbm(spp_list_any, inputs.bodymasses, inputs.is_producer, biomass_adbm, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)

        # 3. Check if ALL models in the group succeeded.
        #    A success is a result that is not `nothing`.
        all_succeeded_bm = (result_ltm !== nothing) && (result_atn !== nothing) && (result_adbm !== nothing)

        # 4. If all 3 succeeded, save them.
        #    If `all_succeeded_bm` is false, we do nothing and loop to the next `run_ID`.
        if all_succeeded_bm
            # We now check each counter *individually* before saving,
            # as one model (e.g., LTM) might finish before another (e.g., ATN).
            
            # Save LTM
            if counters_bm["LTM"] < N_REPLICATES
                counters_bm["LTM"] += 1 # Increment success counter
                fw_id = "ltm_$(counters_bm["LTM"])" # Create success ID
                push!(master_df, (run_ID_bm, fw_id, "LTM", S, target_basal_fraction, missing, result_ltm.percent_basal, result_ltm.connectance, result_ltm.adj, inputs.bodymasses, inputs.metabolic_classes))
            end

            # Save ATN
            if counters_bm["ATN"] < N_REPLICATES
                counters_bm["ATN"] += 1
                fw_id = "atn_$(counters_bm["ATN"])"
                push!(master_df, (run_ID_bm, fw_id, "ATN", S, target_basal_fraction, missing, result_atn.percent_basal, result_atn.connectance, result_atn.adj, inputs.bodymasses, inputs.metabolic_classes))
            end

            # Save ADBM
            if counters_bm["ADBM"] < N_REPLICATES
                counters_bm["ADBM"] += 1
                fw_id = "adbm_$(counters_bm["ADBM"])"
                push!(master_df, (run_ID_bm, fw_id, "ADBM", S, target_basal_fraction, missing, result_adbm.percent_basal, result_adbm.connectance, result_adbm.adj, inputs.bodymasses, inputs.metabolic_classes))
            end
        end
        # --- *** END OF NEW LOGIC *** ---

        # Print a progress update every 500 attempts.
        if run_ID_bm % 500 == 0
            @info "Loop 1 RunID: $run_ID_bm | LTM: $(counters_bm["LTM"]) | ATN: $(counters_bm["ATN"]) | ADBM: $(counters_bm["ADBM"])"
        end
    end
    # --- End of Loop 1 ---
    @info "Loop 1 finished in $run_ID_bm attempts."
    # Warn the user if the loop exited by hitting the max attempt limit.
    if run_ID_bm == MAX_TOTAL_ATTEMPTS_BM
        @warn "Loop 1 hit max attempts ($MAX_TOTAL_ATTEMPTS_BM). Not all replicates may be generated."
    end

    # --- Loop 2: Connectance Models ---
    @info "Running Loop 2 (Connectance Models)..."
    counters_c = Dict("Niche" => 0, "Cascade" => 0, "Random" => 0)
    run_ID_c = 0 # Attempt counter for this loop

    # Same logic as Loop 1.
    while any(values(counters_c) .< N_REPLICATES) && run_ID_c < MAX_TOTAL_ATTEMPTS_C
        run_ID_c += 1
        
        # 1. Generate new shared input for this attempt.
        C_rand = rand(Uniform(CONNECTANCE_RANGE...))

        # --- *** NEW LOGIC: "ALL OR NOTHING" *** ---
        # 2. Run all 3 models in the group FIRST.
        result_niche = run_and_filter_niche(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        result_cascade = run_and_filter_cascade(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        result_random = run_and_filter_random(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)

        # 3. Check if ALL models in the group succeeded.
        all_succeeded_c = (result_niche !== nothing) && (result_cascade !== nothing) && (result_random !== nothing)

        # 4. If (and only if) all 3 succeeded, save them.
        if all_succeeded_c
            # Save Niche
            if counters_c["Niche"] < N_REPLICATES
                counters_c["Niche"] += 1
                fw_id = "niche_$(counters_c["Niche"])"
                push!(master_df, (run_ID_c, fw_id, "Niche", S, missing, C_rand, result_niche.percent_basal, result_niche.connectance, result_niche.adj, missing, missing))
            end
            
            # Save Cascade
            if counters_c["Cascade"] < N_REPLICATES
                counters_c["Cascade"] += 1
                fw_id = "cascade_$(counters_c["Cascade"])"
                push!(master_df, (run_ID_c, fw_id, "Cascade", S, missing, C_rand, result_cascade.percent_basal, result_cascade.connectance, result_cascade.adj, missing, missing))
            end

            # Save Random
            if counters_c["Random"] < N_REPLICATES
                counters_c["Random"] += 1
                fw_id = "random_$(counters_c["Random"])"
                push!(master_df, (run_ID_c, fw_id, "Random", S, missing, C_rand, result_random.percent_basal, result_random.connectance, result_random.adj, missing, missing))
            end
        end
        # --- *** END OF NEW LOGIC *** ---

        # Print progress update.
        if run_ID_c % 500 == 0
            @info "Loop 2 RunID: $run_ID_c | Niche: $(counters_c["Niche"]) | Cascade: $(counters_c["Cascade"]) | Random: $(counters_c["Random"])"
        end
    end
    # --- End of Loop 2 ---
    @info "Loop 2 finished in $run_ID_c attempts."
    if run_ID_c == MAX_TOTAL_ATTEMPTS_C
        @warn "Loop 2 hit max attempts ($MAX_TOTAL_ATTEMPTS_C). Not all replicates may be generated."
    end

else
    # --- 5.B. WORKFLOW A: FIXED ATTEMPTS (For Loops) ---
    # This workflow remains unchanged. It is already designed
    # to run all models for each fixed attempt, which is correct.
    @info "Starting Workflow A: Performing $N_REPLICATES fixed attempts per model. No filtering."
    
    spp_list_int = 1:S
    spp_list_any = convert(Vector{Any}, spp_list_int)

    # Loop exactly `N_REPLICATES` times.
    for i = 1:N_REPLICATES
        if i % (N_REPLICATES/10) == 0
            @info "Running attempt $i / $N_REPLICATES..."
        end

        # --- 1. Generate Inputs for this attempt ---
        target_basal_fraction = rand(Uniform(BASAL_RANGE...))
        inputs = generate_bodymass_inputs(S, target_basal_fraction, BODYMASS_RANGES)
        biomass_adbm = inputs.bodymasses .^ (-0.75)
        C_rand = rand(Uniform(CONNECTANCE_RANGE...))
        
        # In this workflow, run_ID = attempt number, and fw_ID is the "replicate" number.
        run_ID = i 

        # --- 2. Run Body Mass Models ---
        # The wrapper function is called with `verify_web = false`.
        # It will *always* return a result, never `nothing`.
        result_ltm = run_and_filter_ltm(spp_list_int, inputs.bodymasses, inputs.metabolic_classes, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "ltm_$i", "LTM", S, target_basal_fraction, missing, result_ltm.percent_basal, result_ltm.connectance, result_ltm.adj, inputs.bodymasses, inputs.metabolic_classes))
        
        result_atn = run_and_filter_atn(spp_list_any, inputs.bodymasses, inputs.is_producer, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "atn_$i", "ATN", S, target_basal_fraction, missing, result_atn.percent_basal, result_atn.connectance, result_atn.adj, inputs.bodymasses, inputs.metabolic_classes))

        result_adbm = run_and_filter_adbm(spp_list_any, inputs.bodymasses, inputs.is_producer, biomass_adbm, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "adbm_$i", "ADBM", S, target_basal_fraction, missing, result_adbm.percent_basal, result_adbm.connectance, result_adbm.adj, inputs.bodymasses, inputs.metabolic_classes))

        # --- 3. Run Connectance Models ---
        result_niche = run_and_filter_niche(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "niche_$i", "Niche", S, missing, C_rand, result_niche.percent_basal, result_niche.connectance, result_niche.adj, missing, missing))
        
        result_cascade = run_and_filter_cascade(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "cascade_$i", "Cascade", S, missing, C_rand, result_cascade.percent_basal, result_cascade.connectance, result_cascade.adj, missing, missing))

        result_random = run_and_filter_random(S, C_rand, verify_web, BASAL_RANGE, CONNECTANCE_RANGE)
        push!(master_df, (run_ID, "random_$i", "Random", S, missing, C_rand, result_random.percent_basal, result_random.connectance, result_random.adj, missing, missing))
    end
    @info "Workflow A finished."

end # End if verify_web
# --- 6. Save Final Results ---

@info "All replicates generated. Saving results..."

# --- *** Create dynamic filename *** ---
# Create a status string based on the toggle.
status_str = verify_web ? "verified" : "unverified"
# Get the current date.
date_str = Dates.format(now(), "dd-mm-yyyy")
# Combine all parts into the new filename format.
filename_base = "network_test_$(status_str)_seed_$(MASTER_SEED)_$(date_str)"


# Ensure Output Directory Exists
mkpath(OUTPUT_DIR)

# Define Full Paths
csv_path = joinpath(OUTPUT_DIR, "$(filename_base).csv")
jld2_path = joinpath(OUTPUT_DIR, "$(filename_base).jld2")

# --- Save CSV File (Summary) ---
# Create a summary DataFrame that excludes the large objects (matrices, vectors).
# This is for quick inspection.
#@info "Saving summary to CSV..."
#df_summary = select(master_df, Not([:AdjacencyMatrix, :BodyMasses, :MetabolicClasses]))
#try
#    CSV.write(csv_path, df_summary)
#    @info "Successfully saved summary to: $csv_path"
#catch e
#    @warn "Failed to save CSV file!"
#    println(e)
#end

# --- Save CSV File (Full Data) ---
# --- *** CHANGE: Save the full `master_df` to the CSV *** ---
# This will save all columns, including matrices and vectors as strings.
@info "Saving full data to CSV..."
# df_summary = select(master_df, Not([:AdjacencyMatrix, :BodyMasses, :MetabolicClasses])) # This line is removed.
try
    CSV.write(csv_path, master_df) # Changed df_summary to master_df
    @info "Successfully saved full data to: $csv_path"
catch e
    @warn "Failed to save CSV file!"
    println(e)
end

# --- Save JLD2 File (Full Data) ---
# Save the complete `master_df` object, including matrices and vectors,
# to a JLD2 file. This is the file you will use for your analysis.
@info "Saving full data to JLD2..."
try
    # JLD2.save_object is a simple way to save a single variable
    JLD2.save_object(jld2_path, master_df)
    @info "Successfully saved full data (with matrices) to: $jld2_path"
catch e
    @warn "Failed to save JLD2 file!"
    println(e)
end

@info "--- Workflow Complete ---"

