# this files contains functions to process simulation outputs.

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


## -- Function to get alive and connected species --
"""
    get_alive_connected_species(Beq, params)
1. use get_alive_species(Beq) to identify alive species;
2. subset the adjacency matrix to alive species;
3. keep species that have at least one prey or at least one consumer
   among the alive species.
"""
function get_alive_connected_species(Beq, params)
    # Get alive species
    alive_idx = get_alive_species(Beq)

    if isempty(alive_idx)
        return Int[]
    end

    # Original trophic adjacency matrix
    A = params.A

    # Subset to alive species
    A_alive = A[alive_idx, alive_idx]

    # Consumers: rows with at least one prey
    has_prey = vec(sum(A_alive; dims = 2)) .> 0

    # Resources: columns with at least one consumer
    has_consumer = vec(sum(A_alive; dims = 1)) .> 0

    # Keep species that are connected as either consumer or resource
    keep_mask = has_prey .| has_consumer

    alive_connected_idx = alive_idx[keep_mask]

    return alive_connected_idx
end


## -- Function to get the energy flux from prey to predator at equilibrium ----
"""
    energy_flux_classic_formula(B, params; assimilated = false)

Calculate realised trophic fluxes under ClassicResponse using an explicit formula.

For ClassicResponse:

    F[i, j] =
        ω[i, j] * aᵣ[i, j] * B[j]^h /
        (M[i] * (1 + c[i] * B[i] +
         sum_k ω[i, k] * aᵣ[i, k] * hₜ[i, k] * B[k]^h))

Then realised prey-loss flux is:

    flux[i, j] = B[i] * F[i, j]

If assimilated = true, returns:

    flux[i, j] = B[i] * e[i, j] * F[i, j]
"""
function energy_flux_classic_formula(
    B::AbstractVector,
    params;
    assimilated::Bool = false)
    
    fr = params._value.functional_response

    Bsafe = max.(Float64.(B), 0.0)

    ω  = fr.ω
    aᵣ = fr.aᵣ
    hₜ = fr.hₜ
    c  = fr.c
    h  = fr.h
    M  = params.body_mass
    e  = params.e

    S = length(Bsafe)

    I = Int[]
    J = Int[]
    V = Float64[]

    consumers, resources, ωvals = findnz(ω)

    for (i, j, ωij) in zip(consumers, resources, ωvals)

        # Denominator for consumer i
        denom_sum = 0.0
        prey_i = findnz(ω[i, :])[1]

        for k in prey_i
            denom_sum += ω[i, k] * aᵣ[i, k] * hₜ[i, k] * abs(Bsafe[k])^h
        end

        denom = M[i] * (1.0 + c[i] * Bsafe[i] + denom_sum)

        Fij = ωij * aᵣ[i, j] * abs(Bsafe[j])^h / denom

        # prey-loss flux
        flux_ij = Bsafe[i] * Fij

        # predator-gain flux after assimilation
        if assimilated
            flux_ij *= e[i, j]
        end

        push!(I, i)
        push!(J, j)
        push!(V, flux_ij)
    end

    return sparse(I, J, V, S, S)
end

function gini_coefficient(x::AbstractVector)
    x = collect(skipmissing(x))
    x = x[.!isnan.(x)]
    x = x[x .> 0]

    isempty(x) && return NaN
    sum(x) == 0 && return NaN

    x = sort(x)
    n = length(x)

    return 2 * sum((1:n) .* x) / (n * sum(x)) - (n + 1) / n
end


## -- Function to get pairwise interaction strengths at equilibrium --
"""
    get_alive_connected_jacobian(params, Beq, alive_connected_idx)

Return the Jacobian submatrix for alive and connected species.
"""
function get_alive_connected_jacobian(params, Beq, alive_connected_idx)
    isempty(alive_connected_idx) && return Matrix{Float64}(undef, 0, 0)
    J = jacobian(params, Beq)
    return Matrix(J[alive_connected_idx, alive_connected_idx])
end


"""
    get_pairwise_interaction_strengths(params, Beq, alive_connected_idx; absolute = true)

Return off-diagonal Jacobian elements as pairwise interaction strengths.
"""
function get_pairwise_interaction_strengths(
    J::AbstractMatrix;
    absolute::Bool = true
    )
    
    S = size(J, 1)

    S <= 1 && return Float64[]

    offdiag = trues(S, S)
    for i in 1:S
        offdiag[i, i] = false
    end

    vals = J[offdiag]
    vals = vals[.!iszero.(vals)]

    if absolute
        vals = abs.(vals)
    end

    return vals
end

function get_skewness(x::AbstractVector)
    x = collect(skipmissing(x))
    x = x[.!isnan.(x)]

    length(x) < 3 && return NaN

    μ = mean(x)
    σ = std(x)

    σ == 0 && return NaN

    return mean(((x .- μ) ./ σ).^3)
end

## --- Function to calculate resilience -----------------
"""
    get_fw_resilience(params, Beq, alive_connected_idx)

Calculate resilience from the alive-connected Jacobian.
"""
function get_fw_resilience(J::AbstractMatrix)
    size(J, 1) == 0 && return NaN
    return resilience(J)
end

## --- Function to calculate reactivity -----------------
"""
    get_fw_reactivity(params, Beq, alive_connected_idx)

Calculate reactivity from the alive-connected Jacobian.
"""
function get_fw_reactivity(J::AbstractMatrix)
    size(J, 1) == 0 && return NaN
    return reactivity(J)
end


## --- Function to get all output of a simulation -----------------

"""
    get_sim_summary(params, sol)

Return simulation summary metrics at equilibrium.

Returned metrics:
- S_post: number of species in the post-simulation network (for NAN checking)
- L_post: number of links in the post-simulation network (for NAN checking)
- persistence
- biomass_shannon
- gini_fluxes
- skewness_IS
- resilience
- reactivity
- post_adj
"""
function get_sim_summary(params, sol)

    # Original species metabolic classes
    met_class = params.metabolic_class

    # Final state biomasses
    Beq = sol.u[end]

    # Original adjacency matrix
    A = params.A

    # Alive and connected species
    alive_connected_idx = get_alive_connected_species(Beq, params)

    # Alive and connected species metabolic classes
    met_class_alive_connected = met_class[alive_connected_idx]

    # Empty case
    if isempty(alive_connected_idx)
        return (
            S_post = 0,
            L_post = 0,
            persistence = 0.0,
            biomass_shannon = NaN,
            gini_fluxes = NaN,
            skewness_IS = NaN,
            resilience = NaN,
            reactivity = NaN,
            post_adj = spzeros(0, 0),
            met_class_alive_connected = String[]
        )
    end

    # Reduced biomass and adjacency matrix
    Beq_alive_connected = Beq[alive_connected_idx]
    post_adj = A[alive_connected_idx, alive_connected_idx]
    S_post = length(alive_connected_idx)
    L_post = sum(post_adj)

    # Species richness
    S_initial = ENDI.richness(params.A)
    S_eq = length(alive_connected_idx)

    # Persistence
    persistence_eq = S_eq / S_initial

    # Biomass Shannon diversity
    biomass_shannon_eq = ENDI.shannon_diversity(Beq_alive_connected)

    # Fluxes using explicit ClassicResponse formula
    flux_formula = energy_flux_classic_formula(Beq, params; assimilated = false)
    flux_formula_ac = flux_formula[alive_connected_idx, alive_connected_idx]

    flux_formula_vals = nonzeros(sparse(flux_formula_ac))
    flux_formula_vals = flux_formula_vals[flux_formula_vals .> 0]

    gini_fluxes_eq = gini_coefficient(flux_formula_vals)

    # get Jacobian matrix for alive and connected species
    J_alive_connected = get_alive_connected_jacobian(
        params,
        Beq,
        alive_connected_idx
    )

    # get species interaction strengths (off-diagonal Jacobian elements)
    IS_vals = get_pairwise_interaction_strengths(
        J_alive_connected;
        absolute = true
    )

    # Skewness of Jacobian interaction strengths
    skewness_IS_eq = get_skewness(IS_vals)

    # Resilience
    resilience_eq = get_fw_resilience(J_alive_connected)

    # Reactivity
    reactivity_eq = get_fw_reactivity(J_alive_connected)

    return (
        S_post = S_post,
        L_post = L_post,
        persistence = persistence_eq,
        biomass_shannon = biomass_shannon_eq,
        gini_fluxes = gini_fluxes_eq,
        skewness_IS = skewness_IS_eq,
        resilience = resilience_eq,
        reactivity = reactivity_eq,
        post_adj = post_adj,
        met_class_alive_connected = met_class_alive_connected
    )
end

