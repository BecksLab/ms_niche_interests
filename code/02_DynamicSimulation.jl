#=
-------------------------------------------------
02_DynamicSimulation.jl
Runs END to equilibrium on unfiltered networks.
-------------------------------------------------
=#

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using EcologicalNetworksDynamics
using EcoNetPostProcessing
using DifferentialEquations
using Statistics

# --- 2. Load All Code ---
include(joinpath("lib", "sim.jl"));

# --- 3. Import networks .jld2 object ---
networks = load_object("networks/networks.jld2")
sort!(networks, :fw_ID) 

# --- 4. Run Dynamic Simulations ---
tmin, tmax, tspan = 2000, 5000, 500

simulation_summary = DataFrame(
    fw_ID = String[],
    model = String[],
    richness_equilibrium = Float64[],    
    connectance_equilibrium = Float64[], 
    persistence = Float64[],
    max_trophic_level = Float64[],
    total_biomass = Float64[],
    cv_total_biomass = Float64[],
    shannon = Float64[],
    evenness = Float64[],
    resilience = Float64[],
    reactivity = Float64[]
)

final_network = DataFrame(
    fw_ID = String[],
    Model = String[],
    AdjacencyMatrix = Union{Matrix{Int},Missing}[],
    MetabolicClasses = Union{Vector{Symbol},Missing}[] 
)

# --- Execution Loop ---
for i in 1:nrow(networks)
    @info "Simulating food web $(i) / $(nrow(networks))"
    fwid         = networks.fw_ID[i]
    model_name   = networks.Model[i]                
    bodymasses   = networks.BodyMasses[i]
    met_class    = networks.MetabolicClasses[i]
    adj          = networks.AdjacencyMatrix[i]

    fw = Foodweb(adj)  

    if model_name in ("ADBM", "ATN", "LTM")
        bodymasses_rescaled = rescale_bodymass(bodymasses, met_class)
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
            extinction_callback(params, 1e-6; verbose = false),  
            TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass; min_t = tmin),
        ),
        show_degenerated = false,
    )

    # Process and safely trap potential numerical errors
    out = try
        get_sim_summary(params, sol, tspan)
    catch err
        @warn "Simulation failed for fw_ID=$fwid" err
        (
            richness_equilibrium = NaN,
            connectance_equilibrium = NaN,
            persistence = NaN,
            max_trophic_level = NaN,
            total_biomass = NaN,
            cv_total_biomass = NaN,
            shannon = NaN,
            evenness = NaN,
            resilience = NaN,
            reactivity = NaN,
            alive_connected_A = missing
        )
    end

    # Package output objects for storage
    summary_out = NamedTuple{filter(k -> k != :alive_connected_A, keys(out))}(out)
    push!(simulation_summary, (; fw_ID=string(fwid), model=string(model_name), summary_out...); promote = true)
    push!(final_network, (fw_ID=string(fwid), Model=string(model_name), AdjacencyMatrix=Int.(out.alive_connected_A), MetabolicClasses=met_class))
end

# --- 5. Save Simulation Summary ---
CSV.write("outputs/summary_stability.csv", simulation_summary)
JLD2.save_object("networks/networks_END.jld2", final_network)