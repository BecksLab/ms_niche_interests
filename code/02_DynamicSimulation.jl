#=
This script performs dynamic simulations of generated food webs 
and generates two output files.

First, a jld2 file (networks_END.jld2) containing the post-simulation food web networks, 
       represented as adjacency matrices of feeding links among species that remain alive and connected at equilibrium.

Second, a CSV file (dynamic_metrics.csv) containing dynamic stability-related metrics for each food web:
        - S_post: number of species in the post-simulation network (save for NAN checking)
        - L_post: number of links in the post-simulation network (save for NAN checking)
        - Persistence; 
        - Shannon diversity of biomass distribution; 
        - Gini coefficient of energy fluxes;
        - Skewness of absolute Jacobian interaction strengths;
        - Resilience;
        - Reactivity; 
=#

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using EcologicalNetworksDynamics
using EcoNetPostProcessing
using DifferentialEquations
using DiffEqCallbacks
using Statistics
using SparseArrays

# --- 2. Load All Code ---
include(joinpath("lib", "sim.jl"));

# --- 3. Import networks .jld2 object ---
pre_networks = load_object("networks/networks.jld2")
sort!(pre_networks, :fw_ID) 

# --- 4. Run Dynamic Simulations ---
tmin, tmax = 2000, 5000

dynamic_metrics = DataFrame(
    fw_ID = String[],
    Model = String[],
    S_post = Int64[], # Number of species in the post-simulation network
    L_post = Int64[], # Number of links in the post-simulation network
    persistence = Float64[], # Persistence
    biomass_shannon = Float64[], # Shannon diversity of biomass distribution
    gini_fluxes = Float64[], # Gini coefficient of energy fluxes
    skewness_IS = Float64[], # Skewness of absolute interaction strengths
    resilience = Float64[],  # Resilience
    reactivity = Float64[]  # Reactivity
)

post_networks = DataFrame(
    fw_ID = String[],
    Model = String[],
    AdjacencyMatrix = Any[]
)

for i in 1:nrow(pre_networks)
    @info "Simulating food web $(i) / $(nrow(pre_networks))"
    fwid         = pre_networks.fw_ID[i]
    model_name   = pre_networks.Model[i]                
    body_masses  = pre_networks.BodyMasses[i]
    met_class    = pre_networks.MetabolicClasses[i]
    pre_adj      = pre_networks.AdjacencyMatrix[i]


    fw = Foodweb(pre_adj)  

    if model_name in ("ADBM", "ATN", "LTM")
        bodymasses_rescaled = rescale_bodymass(body_masses, met_class)
        params = default_model(
            fw,
            BodyMass(bodymasses_rescaled),
            ClassicResponse(; h = 2),
        )
    else
        # Topology-only models revert to standard uniform mass assumptions
        params = default_model(
            fw,
            BodyMass(; Z = 100),
            ClassicResponse(; h = 2),
        )
    end

    # Set initial uniform random biomasses
    B0 = rand(params.S)

    # Execute simulation with performance callbacks active
    sol = simulate(
        params, B0, tmax;
        callback = CallbackSet(
            extinction_callback(params, 1e-6; verbose = false),  # extinction threshold 1e-6
            TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass; min_t = tmin),
        ),
        show_degenerated = false,
    )

    # Process and safely trap potential numerical errors
    out = try
        get_sim_summary(params, sol)
    catch err
        @warn "Simulation failed for fw_ID=$fwid" err
        (
            persistence = NaN,
            S_post = NaN,
            L_post = NaN,
            biomass_shannon = NaN,
            gini_fluxes = NaN,
            skewness_IS = NaN,
            resilience = NaN,
            reactivity = NaN,
            post_adj = missing
        )
    end

    # Remove matrix object before saving summary to CSV
    summary_out = NamedTuple{
        filter(k -> k != :post_adj, keys(out))
    }(out)

    # Save scalar dynamic metrics.
    push!(dynamic_metrics, (; fw_ID=string(fwid), Model = string(model_name), summary_out...); 
                            promote = true)

    # Save the post-simulation adjacency matrix.
    push!(post_networks, (; fw_ID = string(fwid), Model = string(model_name), AdjacencyMatrix = Int.(out.post_adj));
                          promote = true)
end


# --- 5. Save Simulation Summary ---
CSV.write("outputs/dynamic_metrics.csv", dynamic_metrics)
JLD2.save_object("networks/networks_END.jld2", post_networks)
