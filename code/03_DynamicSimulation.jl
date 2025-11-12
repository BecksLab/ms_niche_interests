#=
Dynamic simulations are conducted in END.
Simulation output include:
- Persistence; 
- Maximum trophic level; 
- Total Biomass; 
- CV of total biomass;  
- Shannon; 
- Evenness; 
- Resilience if at steady state; (using EcoNetPostProcessing.jl) 
- Reactivity; (using EcoNetPostProcessing.jl) 
- Robustness; (using EcoNetPostProcessing.jl)
=#

using Pkg
Pkg.activate(".")

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using EcologicalNetworksDynamics
using EcoNetPostProcessing
using DifferentialEquations
using Statistics

# --- 2. Load All Code ---
include("lib/sim.jl");

# --- 3. Import networks .jld2 object ---
path = joinpath(@__DIR__, "data", "outputs", "network_test_verified_seed_42_29-10-2025.jld2")
networks = load_object(path)
sort!(networks, :fw_ID) 

# --- 4. Run Dynamic Simulations ---
tmin, tmax, tspan = 3000, 6000, 1000

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

for i in 1:10 # nrow(networks)
    fwid         = networks.fw_ID[i]
    model_name   = networks.Model[i]                
    bodymasses   = networks.BodyMasses[i]
    met_class    = networks.MetabolicClasses[i]
    adj            = networks.AdjacencyMatrix[i]

    fw = Foodweb(adj)  

    if model_name in ("ADBM", "ATN", "LTM")
        bodymasses_rescaled = rescale_bodymass(bodymasses, met_class)
        params = default_model(
            fw,
            BodyMass(bodymasses_rescaled),
            ClassicResponse(; h = 2),
        )
    else
        params = default_model(
            fw,
            BodyMass(; Z = 100),
            ClassicResponse(; h = 2),
        )
    end

    # Choose an initial biomass (adjust to your API if you have a helper for this)
    B0 = rand(params.S) 

    sol = simulate(
        params, B0, tmax;
        callback = CallbackSet(
            extinction_callback(params, 1e-6; verbose = false),  # extinction threshold
            TerminateSteadyState(1e-8, 1e-6, DiffEqCallbacks.allDerivPass; min_t = tmin),
        ),
        show_degenerated = false,
    )

    out = get_sim_summary(params, sol, tspan)

    push!(simulation_summary, (
        fw_ID = fwid,
        model = model_name,
        out...
    ))
end

# --- 6. Save Simulation Summary ---
CSV.write("data/outputs/END_simulation.csv", simulation_summary)