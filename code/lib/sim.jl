### This file contains functions to be used in the pre-perturbation simulation.
ENDI = EcologicalNetworksDynamics.Internals

function rescale_bodymass(bodymasses::AbstractVector, met_class::AbstractVector)
    idx = findall(==(:producer), met_class)
    ref = minimum(bodymasses[idx])
    bodymasses ./= float(ref)
    return bodymasses
end

function get_spe_generality(A)
    A = sparse(A)
    S = size(A, 1)
    generality = zeros(S)
    for i in 1:S
      preys = ENDI.preys_of(i, A)
      length_preys = length(preys)
      generality[i] = length_preys
    end
    return generality
end

function get_spe_vulnerability(A)
    A = sparse(A)
    S = size(A, 1)
    vulnerability = zeros(S)
    for i in 1:S
      predators = ENDI.predators_of(i, A)
      length_predators = length(predators)
      vulnerability[i] = length_predators
    end
    return vulnerability
end

function get_fw_generality(A::AbstractMatrix)
    g = get_spe_generality(A)
    nz = g[g .> 0]
    if isempty(nz)
        return NaN
    end
    mean_generality = sum(g) / length(nz)   # same as R: sum(colSums)/count(colSum>0)
    sd_generality   = length(nz) > 1 ? std(nz) : 0.0
    return mean_generality, sd_generality
end

function get_fw_vulnerability(A::AbstractMatrix)
    g = get_spe_vulnerability(A)
    nz = g[g .> 0]
    if isempty(nz)
        return NaN
    end
    mean_vulnerability = sum(g) / length(nz)   # same as R: sum(colSums)/count(colSum>0)
    sd_vulnerability   = length(nz) > 1 ? std(nz) : 0.0
    return mean_vulnerability, sd_vulnerability
end

function get_IS_summary(params, Beq, alive_idx)
    j = jacobian(params, Beq)
    j_alive = j[alive_idx, alive_idx]
    j_alive_abs = abs.(j_alive)
    # all top-down interaction strengths
    mean_abs_IS_all = mean(j_alive_abs[j_alive_abs .> 0])
    sd_abs_IS_all = std(j_alive_abs[j_alive_abs .> 0])
    skew_abs_IS_all = skewness(j_alive_abs[j_alive_abs .> 0])  
    top_down_IS = tril(j_alive_abs, -1)
    mean_top_down = mean(top_down_IS[top_down_IS .> 0])
    sd_top_down = std(top_down_IS[top_down_IS .> 0])
    bottom_up_IS = triu(j_alive_abs, 1)
    mean_bottom_up = mean(bottom_up_IS[bottom_up_IS .> 0])
    sd_bottom_up = std(bottom_up_IS[bottom_up_IS .> 0])
    IS_summary = (
         mean_top_down = mean_top_down,
         mean_bottom_up = mean_bottom_up,
         mean_abs_IS_all = mean_abs_IS_all,
         sd_top_down = sd_top_down,
         sd_bottom_up = sd_bottom_up,
         sd_abs_IS_all = sd_abs_IS_all,
         skew_abs_IS_all = skew_abs_IS_all
      )
    return IS_summary
end

function get_fw_resilience(params, Beq, alive_idx)
    j = jacobian(params, Beq)
    j_alive = j[alive_idx, alive_idx]
    lumbda_eq = resilience(j_alive)
    return lumbda_eq
end

function get_fw_reactivity(params, Beq, alive_idx)
    j = jacobian(params, Beq)
    j_alive = j[alive_idx, alive_idx]
    reactivity_eq = reactivity(j_alive)
    return reactivity_eq
end


## --- Function to get all output of a simulation -----------------
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
    richness_equilibrium = S_eq,
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

