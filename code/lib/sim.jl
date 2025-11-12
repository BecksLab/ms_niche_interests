### This file contains functions to be used in the pre-perturbation simulation.
ENDI = EcologicalNetworksDynamics.Internals

## --- Function to rescale bodymass -----------------
"""
    rescale_bodymass(bodymasses::AbstractVector, met_class::AbstractVector)
    bodymasses = vector of bodymasses
    met_class = vector of metabolic classes (:producer, :herbivore, :carnivore, etc.)  
    Rescale bodymasses so that the smallest producer has bodymass = 1. 
"""
function rescale_bodymass(bodymasses::AbstractVector, met_class::AbstractVector)
    idx = findall(==(:producer), met_class)
    ref = minimum(bodymasses[idx])
    bodymasses ./= float(ref)
    return bodymasses
end


## --- Function to calculate resilience -----------------
"""
    get_fw_resilience(params::ModelParameters, Beq::AbstractVector, alive_idx::AbstractVector)
    params = model parameters
    Beq = equilibrium biomasses
    alive_idx = indices of alive species
    Calculate the resilience of the food web at equilibrium.
"""
function get_fw_resilience(params, Beq, alive_idx)
    j = jacobian(params, Beq)
    j_alive = j[alive_idx, alive_idx]
    lumbda_eq = resilience(j_alive)
    return lumbda_eq
end

## --- Function to calculate reactivity -----------------
"""
    get_fw_reactivity(params::ModelParameters, Beq::AbstractVector, alive_idx::AbstractVector)
    params = model parameters
    Beq = equilibrium biomasses
    alive_idx = indices of alive species
    Calculate the reactivity of the food web at equilibrium.
""" 
function get_fw_reactivity(params, Beq, alive_idx)
    j = jacobian(params, Beq)
    j_alive = j[alive_idx, alive_idx]
    reactivity_eq = reactivity(j_alive)
    return reactivity_eq
end


## --- Function to get all output of a simulation -----------------
"""
    get_sim_summary(params::ModelParameters, sol::ODESolution, tspan::Int64)
    params = model parameters
    sol = solution of the ODE simulation
    tspan = number of timesteps to consider at the end of the simulation for time series metrics
    Return a named tuple with various summary statistics of the simulation at equilibrium.
"""
function get_sim_summary(params, sol, tspan)
  # Get bomass dynamics equally distributed at last tspan timesteps using interpolation
   last_t = sol.t[end]
   vals = sol(last_t - tspan + 1:1:last_t).u
   biomass_mat = transpose(hcat(vals...)) # tspan x species

   # final state biomasses
   Beq = sol.u[end]

   # --- Get alive and connected species ---
   alive_idx = get_alive_species(Beq)
   A = params.A
   # find alive producers (rows all zeros = no prey) and consumers (cols all zeros = no consumers)
   has_prey = vec(sum(A[alive_idx, alive_idx]; dims=2)) .> 0
   has_consumer = vec(sum(A[alive_idx, alive_idx]; dims=1)) .> 0
   # keep species that are alive and have at least one prey or one consumer
   keep_mask = has_prey .| has_consumer
   alive_idx = alive_idx[keep_mask]
   # reduced vectors/matrices
   Beq_alive = Beq[alive_idx]
   alive_A   = A[alive_idx, alive_idx]

   # --- for cases where all species go extinct or become disconnected ---
    if isempty(alive_idx)
        return (
            richness_equilibrium = 0,
            connectance_equilibrium = 0.0,
            persistence = 0.0,
            max_trophic_level = 0.0,
            total_biomass = 0.0,
            cv_total_biomass = 0.0,
            shannon = 0.0,
            evenness = 0.0,
            resilience = NaN,
            reactivity = NaN
        )
    end

   # species richness and connectance at equilibrium
   S_eq = length(alive_idx)
   C_eq = round(count(!iszero, alive_A) / (S_eq * S_eq), digits=3)

   # maximum trophic level at pre-perturbation equilibrium
   tls = ENDI.trophic_levels(alive_A)
   maxTL_eq = maximum(tls)  

   # total biomass at equilibrium
   total_biomass_eq = sum(Beq_alive)
   
   # CV of total biomass
   B_mean_eq = mean(sum(biomass_mat[:, alive_idx]; dims=2))
   B_sd_eq   = std(sum(biomass_mat[:, alive_idx]; dims=2))
   CV_eq    = B_sd_eq / B_mean_eq
   
   # shannon diversity of biomass distribution
   biomass_shannon_eq = ENDI.shannon_diversity(Beq_alive)

   # evenness of biomass distribution
   biomass_evenness_eq = ENDI.evenness(Beq_alive)

   # resilience
   lumbda_eq = get_fw_resilience(params, Beq, alive_idx)

   # reactivity
    reactivity_eq = get_fw_reactivity(params, Beq, alive_idx)

  return(
    richness_equilibrium = Int(S_eq),
    connectance_equilibrium = C_eq,
    persistence = S_eq / params.S,
    max_trophic_level = maxTL_eq,
    total_biomass = total_biomass_eq,
    cv_total_biomass = CV_eq,
    shannon = biomass_shannon_eq,
    evenness = biomass_evenness_eq,
    resilience = lumbda_eq,
    reactivity = reactivity_eq
    )
end

