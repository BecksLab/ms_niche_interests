#=
-------------------------------------------------
lib/verifying_networks.jl
Helper library for the generative model comparison project.
-------------------------------------------------

This file contains all functional logic for the main build script:
1.  Input Generators (e.g., `generate_bodymass_inputs`)
2.  Property Calculators (e.g., `calculate_connectance`)
3.  Property Checkers (e.g., `is_in_basal_range`)
4.  Internal Processing (e.g., `_process_web`)
5.  Model Wrappers (e.g., `run_and_filter_ltm`)
=#

# --- 1. Dependencies ---
# These are assumed to be loaded by the main script
using Distributions, Random, Statistics, LinearAlgebra

# --- 2. Input Generator ---

"""
    generate_bodymass_inputs(S, target_basal_fraction, bm_ranges)

Generates a consistent set of body mass, metabolic class, and initial biomass
vectors for the body-mass-based models.

Body mass ranges are provided as log10 exponents.

# Arguments
- `S::Int`: Species richness.
- `target_basal_fraction::Float64`: The target fraction of producers (e.g., 0.2).
- `bm_ranges::NamedTuple`: A named tuple with producer/invertebrate log10 ranges.

# Returns
- `NamedTuple`: Contains `bodymasses`, `metabolic_classes`, and `is_producer`.
"""
function generate_bodymass_inputs(
    S::Int,
    target_basal_fraction::Float64,
    bm_ranges::NamedTuple,
)
    # --- 1. Determine Metabolic Classes ---
    # Round the number of producers based on the target fraction.
    num_producers = round(Int, S * target_basal_fraction)
    # Ensure at least one producer if target > 0, and not more than S.
    if num_producers == 0 && target_basal_fraction > 0
        num_producers = 1
    elseif num_producers > S
        num_producers = S
    end
    num_consumers = S - num_producers

    # Create the vector of metabolic class labels.
    metabolic_classes = vcat(
        fill(:producer, num_producers),
        fill(:invertebrate, num_consumers),
    )
    # Randomly shuffle the assignments.
    shuffle!(metabolic_classes)
    
    # --- 2. Create Boolean Vector ---
    # Create the boolean vector needed by some models.
    # We must explicitly convert the BitVector (from .==) to a Vector{Bool}
    # to ensure type compatibility with the lmatrix function.
    is_producer::Vector{Bool} = (metabolic_classes .== :producer)

    # --- 3. Generate Body Masses ---
    bodymasses = zeros(Float64, S)
    # Define log-uniform distributions from the log10 ranges.
    # 10^bm_ranges.producer.min converts the exponent to an absolute value.
    producer_dist = LogUniform(10^bm_ranges.producer.min, 10^bm_ranges.producer.max)
    consumer_dist = LogUniform(10^bm_ranges.invertebrate.min, 10^bm_ranges.invertebrate.max)

    # Iterate and assign a random body mass from the correct distribution.
    for i = 1:S
        if is_producer[i]
            bodymasses[i] = rand(producer_dist)
        else
            bodymasses[i] = rand(consumer_dist)
        end
    end

    # Return all generated inputs.
    return (
        bodymasses = bodymasses,
        metabolic_classes = metabolic_classes,
        is_producer = is_producer,
    )
end


# --- 3. Property Calculators ---

"""
    calculate_connectance(adj::AbstractMatrix) -> Float64

Calculates the directed connectance (C) of a food web.
C = L / S^2.

# Arguments
- `adj::AbstractMatrix`: An adjacency matrix (S x S).

# Returns
- `Float64`: The connectance value.
"""
function calculate_connectance(adj::AbstractMatrix)
    S = size(adj, 1)
    if S == 0
        return 0.0
    end
    # Sum all elements in the matrix to get the total number of links (L).
    L = sum(adj)
    # Divide by the total possible links, S*S.
    return L / (S^2)
end

"""
    calculate_emergent_producers(adj::AbstractMatrix) -> Float64

Calculates the fraction of species that are basal (producers).
A basal species has an out-degree of 0 (e.g., its row sums to 0).

# Arguments
- `adj::AbstractMatrix`: An adjacency matrix (S x S).

# Returns
- `Float64`: The fraction of species that are basal (e.g., 0.2).
"""
function calculate_emergent_producers(adj::AbstractMatrix)
    S = size(adj, 1)
    if S == 0
        return 0.0
    end
    # Sum along rows (dims=2) to get the out-degree for each species.
    out_degrees = sum(adj, dims = 2) # S x 1 Matrix
    # Count how many species have an out-degree of 0.
    basal_count = sum(vec(out_degrees) .== 0)
    # Return the fraction.
    return basal_count / S
end


# --- 4. Property Checkers ---

"""
    is_in_basal_range(percent_basal, basal_range) -> Bool

Checks if the emergent basal percentage is within the allowed range.

# Arguments
- `percent_basal::Float64`: The emergent basal fraction.
- `basal_range::Tuple{Float64,Float64}`: (min, max) allowed range.

# Returns
- `Bool`: `true` if in range, `false` otherwise.
"""
function is_in_basal_range(percent_basal::Float64, basal_range::Tuple{Float64,Float64})
    # Check if value is between min (index 1) and max (index 2).
    return basal_range[1] <= percent_basal <= basal_range[2]
end

"""
    is_in_connectance_range(connectance, connectance_range) -> Bool

Checks if the emergent connectance is within the allowed range.

# Arguments
- `connectance::Float64`: The emergent connectance.
- `connectance_range::Tuple{Float64,Float64}`: (min, max) allowed range.

# Returns
- `Bool`: `true` if in range, `false` otherwise.
"""
function is_in_connectance_range(connectance::Float64, connectance_range::Tuple{Float64,Float64})
    # Check if value is between min (index 1) and max (index 2).
    return connectance_range[1] <= connectance <= connectance_range[2]
end


# --- 5. Internal Processing Function ---

"""
    _process_web(adj, verify_web, basal_range, connectance_range)

Internal helper. Calculates properties of a generated web and
applies filters based on the `verify_web` toggle.

# Arguments
- `adj::Matrix{Int}`: The generated adjacency matrix.
- `verify_web::Bool`: Toggle to enable/disable filtering.
- `basal_range::Tuple{Float64,Float64}`: (min, max) allowed basal %.
- `connectance_range::Tuple{Float64,Float64}`: (min, max) allowed connectance.

# Returns
- `NamedTuple`: (adj, percent_basal, connectance) on success or if `verify_web` is false.
- `nothing`: If `verify_web` is true and the web fails a filter.
"""
function _process_web(
    adj::Matrix{Int},
    verify_web::Bool,
    basal_range::Tuple{Float64,Float64},
    connectance_range::Tuple{Float64,Float64},
)
    # --- 1. Calculate Emergent Properties ---
    # Always calculate properties, regardless of verification,
    # so we can save them in Workflow A.
    percent_basal = calculate_emergent_producers(adj)
    connectance = calculate_connectance(adj)

    # --- 2. Apply Filters (if verify_web is true) ---
    if verify_web
        # --- *** NEW LOGIC (User Request) *** ---
        # All models are checked against BOTH filters.

        # Check 1: Emergent Basal %
        if !is_in_basal_range(percent_basal, basal_range)
            return nothing # Failed verification
        end
        
        # Check 2: Emergent Connectance
        if !is_in_connectance_range(connectance, connectance_range)
            return nothing # Failed verification
        end
        # --- *** END OF NEW LOGIC *** ---
    end

    # --- 3. Return Success ---
    # If verify_web is false, or if it's true and both checks passed,
    # return the results.
    return (adj = adj, percent_basal = percent_basal, connectance = connectance)
end


# --- 6. Wrapper Functions ---
# These functions standardise the interface for all 6 models.
# They call the core model, get the `adj` matrix, and pass it
# to the `_process_web` function for filtering and calculation.

"""
    run_and_filter_ltm(...)

Wrapper for the Latent Trait Model.
"""
function run_and_filter_ltm(
    species_indices::AbstractVector{Int},
    bodymasses::Vector{Float64},
    metabolic_classes::AbstractVector{Symbol},
    verify_web::Bool,
    basal_range::Tuple{Float64,Float64},
    connectance_range::Tuple{Float64,Float64},
)
    # 1. Call the core model function.
    ltm_result = ltm(species_indices, bodymasses, metabolic_classes)
    # 2. Process and filter the resulting matrix.
    #    (Note: model_type argument removed from _process_web)
    return _process_web(
        ltm_result.binary_matrix,
        verify_web,
        basal_range,
        connectance_range,
    )
end

"""
    run_and_filter_atn(...)

Wrapper for the ATN (L-Matrix) Model.
"""
function run_and_filter_atn(
    spp_list_any::Vector{Any},
    bodymass::Vector{Float64},
    is_producer::Vector{Bool},
    verify_web::Bool,
    basal_range::Tuple{Float64,Float64},
    connectance_range::Tuple{Float64,Float64},
)
    # 1. Call the core model function (converts Bool matrix to Int).
    adj = Int.(lmatrix(spp_list_any, bodymass, is_producer))
    # 2. Process and filter the resulting matrix.
    return _process_web(adj, verify_web, basal_range, connectance_range)
end

"""
    run_and_filter_adbm(...)

Wrapper for the Allometric Diet Breadth Model.
"""
function run_and_filter_adbm(
    spp_list_any::Vector{Any},
    bodymass::Vector{Float64},
    is_producer::Vector{Bool},
    biomass_adbm::Vector{Float64},
    verify_web::Bool,
    basal_range::Tuple{Float64,Float64},
    connectance_range::Tuple{Float64,Float64},
)
    # 1. Get model-specific parameters.
    params = adbm_parameters(bodymass, is_producer; Nmethod = :biomass)
    # 2. Call the core model function (converts Bool matrix to Int).
    adj = Int.(adbm(spp_list_any, params, biomass_adbm))
    # 3. Process and filter the resulting matrix.
    return _process_web(adj, verify_web, basal_range, connectance_range)
end

"""
    run_and_filter_niche(...)

Wrapper for the Niche Model.
"""
function run_and_filter_niche(
    S::Int,
    C_rand::Float64,
    verify_web::Bool,
    basal_range::Tuple{Float64,Float64},
    connectance_range::Tuple{Float64,Float64},
)
    # 1. Call the core model function.
    #    This assumes `generate_niche_model` is the *modified* version
    #    that does no internal filtering and just returns an adj matrix.
    adj = generate_niche_model(S, C_rand)
    # 2. Process and filter the resulting matrix.
    return _process_web(adj, verify_web, basal_range, connectance_range)
end

"""
    run_and_filter_cascade(...)

Wrapper for the Cascade Model.
"""
function run_and_filter_cascade(
    S::Int,
    C_rand::Float64,
    verify_web::Bool,
    basal_range::Tuple{Float64,Float64},
    connectance_range::Tuple{Float64,Float64},
)
    # 1. Call the core model function.
    #    This assumes `generate_cascade_model` is the *modified* version.
    adj = generate_cascade_model(S, C_rand)
    # 2. Process and filter the resulting matrix.
    return _process_web(adj, verify_web, basal_range, connectance_range)
end

"""
    run_and_filter_random(...)

Wrapper for the Random Model.
"""
function run_and_filter_random(
    S::Int,
    C_rand::Float64,
    verify_web::Bool,
    basal_range::Tuple{Float64,Float64},
    connectance_range::Tuple{Float64,Float64},
)
    # 1. Call the core model function.
    #    This assumes `generate_random_model` is the *modified* version.
    adj = generate_random_model(S, C_rand)
    # 2. Process and filter the resulting matrix.
    return _process_web(adj, verify_web, basal_range, connectance_range)
end

