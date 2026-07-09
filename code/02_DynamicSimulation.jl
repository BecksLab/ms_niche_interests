#=
This script performs dynamic simulations of generated food webs 
and generates two output files.

First, a jld2 file (networks_END.jld2) containing：
    - burn-in final network structure: AdjacencyMatrix_burn
    - final equilibrium network structure: AdjacencyMatrix
    - final metabolic classes

Second, a CSV file (stability_metrics.csv) containing dynamic stability-related metrics:
    - S_post: number of species in the post-simulation network (save for NAN checking)
    - L_post: number of links in the post-simulation network (save for NAN checking)
    - Persistence: ratio of richness in the final post-simulation network to the original pre-simulation network
    - Shannon diversity of biomass distribution
    - Gini coefficient of consumption fluxes
    - Skewness of absolute off-diagonal Jacobian interaction strengths
    - Resilience
    - Reactivity
    - run_issues: any issues encountered during the simulation (e.g. solver errors, non-steady state, etc.)

Simulation setup:
1. Burn-in run: simulate the full pre-network. Note that some species may remain alive but 
    slowly decline, e.g. consumers with no remaining prey.
2. Filtering: retain only alive species with biologically meaningful links: 
    - initial producers must still have at least one alive consumer, and
    - initial consumers must still have at least one alive prey.
3. Equilibrium run: re-simulate the filtered network to steady state. 
    The final biomass and network are then used to calculate the Jacobian 
    and dynamic stability metrics.
    
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
using Random

Random.seed!(2026)

# --- 2. Load All Code ---
include(joinpath("lib", "sim.jl"));

# --- 3. Import pre-networks .jld2 object ---
pre_networks = load_object("networks/networks.jld2")
sort!(pre_networks, :fw_ID)

# --- 4. Run Dynamic Simulations ---
tmin, tmax = 2000, 20000

stability_metrics = DataFrame(
    fw_ID=String[],
    Model=String[],
    S_post=Int64[], # Number of species in the post-simulation network
    L_post=Int64[], # Number of links in the post-simulation network
    persistence=Float64[], # Persistence
    biomass_shannon=Float64[], # Shannon diversity of biomass distribution
    gini_fluxes=Float64[], # Gini coefficient of consumption fluxes
    skewness_IS=Float64[], # Skewness of absolute off-diagonal Jacobian interaction strengths
    resilience=Float64[],  # Resilience
    reactivity=Float64[],  # Reactivity
    run_issues=String[]
)

post_networks = DataFrame(
    fw_ID=String[],
    Model=String[],
    AdjacencyMatrix_burn=Any[],
    AdjacencyMatrix=Any[],
    MetabolicClasses=Any[]
)

function failed_out(issue; burn_adj = missing)
    return (
        persistence = missing,
        S_post = missing,
        L_post = missing,
        biomass_shannon = missing,
        gini_fluxes = missing,
        skewness_IS = missing,
        resilience = missing,
        reactivity = missing,
        post_adj = missing,
        met_class_post = missing,
        burn_adj = burn_adj,
        run_issues = issue
    )
end

for i in 1:nrow(pre_networks)
    println("Simulating food web $(i) / $(nrow(pre_networks))")

    fwid = pre_networks.fw_ID[i]
    model_name = pre_networks.Model[i]
    body_masses = pre_networks.BodyMasses[i]
    met_class = pre_networks.MetabolicClasses[i]
    pre_adj = pre_networks.AdjacencyMatrix[i]
    S_ori = pre_networks.S[i]

    fw = Foodweb(pre_adj)

    if model_name in ("ADBM", "ATN", "LTM")
        bodymasses_rescaled = rescale_bodymass(body_masses, met_class)
        params_burn = default_model(
            fw,
            BodyMass(bodymasses_rescaled),
            ClassicResponse(; h=2),
        )
    else
        # Topology-only models revert to standard uniform mass assumptions
        params_burn = default_model(
            fw,
            BodyMass(; Z=100),
            ClassicResponse(; h=2),
        )
    end

    # Set initial uniform random biomasses
    B0 = rand(params_burn.S)

    # Run the burn-in simulation to filter out species that are not biologically meaningful
    # Note that some cases already got error in sol_burn due to solver numerical issues
    sol_burn_final_state = try
        sol_burn = simulate(
            params_burn, B0, tmax;
            callback=CallbackSet(
                extinction_callback(params_burn, 1e-6; verbose=false),
                TerminateSteadyState(1e-14, 1e-12, DiffEqCallbacks.allDerivPass, min_t=tmin),
            ),
            show_degenerated=false
        )

        (
            summary=get_sim_summary(params_burn, sol_burn, S_ori),
            issue=""
        )

    catch err
        msg = sprint(showerror, err)
        println(" ==> Burn-in simulation failed for fw_ID=$fwid")
        println("   error = ", msg)

        (
            summary=missing,
            issue="burn_solver_failed: $msg"
        )
    end

    out = if ismissing(sol_burn_final_state.summary)

        failed_out(sol_burn_final_state.issue)

        # If burn-in already collapsed the network, skip the second simulation.
    elseif sol_burn_final_state.summary.S_post == 0 || sol_burn_final_state.summary.L_post == 0
        burn_adj = ismissing(sol_burn_final_state.summary.post_adj) ?
                   missing :
                   Int.(sol_burn_final_state.summary.post_adj)

        (
            persistence=sol_burn_final_state.summary.persistence,
            S_post=sol_burn_final_state.summary.S_post,
            L_post=sol_burn_final_state.summary.L_post,
            biomass_shannon=missing,
            gini_fluxes=missing,
            skewness_IS=missing,
            resilience=missing,
            reactivity=missing,
            post_adj=missing,
            met_class_post=missing,
            burn_adj=burn_adj,
            run_issues="burn_collapsed"
        )

    else
        # Save the burn-in filtered network structure.
        burn_adj = Int.(sol_burn_final_state.summary.post_adj)

        try
            params = default_model(
                Foodweb(sol_burn_final_state.summary.post_adj),
                BodyMass(sol_burn_final_state.summary.body_mass_post),
                MetabolicClass(sol_burn_final_state.summary.met_class_post),
                ClassicResponse(; h=2)
            )

            sol = simulate(
                params, sol_burn_final_state.summary.Beq_post, tmax;
                callback=CallbackSet(
                    extinction_callback(params, 1e-6; verbose=false),
                    TerminateSteadyState(1e-14, 1e-12, DiffEqCallbacks.allDerivPass, min_t=tmin),
                ),
                show_degenerated=false
            )

            eq = check_steady_state(
                sol;
                n_last=5, dt=1.0, extinction_threshold=1e-6, cvtol=1e-3 )

            if !eq.steady
                println(" ==> Food web did not reach biological steady state fw_ID=$fwid")
                println("   max_cv = ", eq.max_cv)
                failed_out(
                    "final_not_steady; maxCV=$(eq.max_cv)";
                    burn_adj=burn_adj
                )
            else
                final_summary = get_sim_summary(params, sol, S_ori)
                (
                    persistence=final_summary.persistence,
                    S_post=final_summary.S_post,
                    L_post=final_summary.L_post,
                    biomass_shannon=final_summary.biomass_shannon,
                    gini_fluxes=final_summary.gini_fluxes,
                    skewness_IS=final_summary.skewness_IS,
                    resilience=final_summary.resilience,
                    reactivity=final_summary.reactivity,
                    post_adj=final_summary.post_adj,
                    met_class_post=final_summary.met_class_post,
                    burn_adj=burn_adj,
                    run_issues="success"
                )
            end

        catch err
            msg = sprint(showerror, err)
            println(" ==> Final simulation failed for fw_ID=$fwid")
            println("   error = ", msg)

            failed_out("final_solver_failed: $msg"; burn_adj=burn_adj)
        end
    end

    # Save scalar dynamic metrics.
    push!(stability_metrics, (
            fw_ID=string(fwid),
            Model=string(model_name),
            S_post=out.S_post,
            L_post=out.L_post,
            persistence=out.persistence,
            biomass_shannon=out.biomass_shannon,
            gini_fluxes=out.gini_fluxes,
            skewness_IS=out.skewness_IS,
            resilience=out.resilience,
            reactivity=out.reactivity,
            run_issues = out.run_issues
        ); promote=true)

    # Save burn_in and final adjacency matrix.
    burn_adj_out = ismissing(out.burn_adj) ? missing : Int.(out.burn_adj)
    post_adj_out = ismissing(out.post_adj) ? missing : Int.(out.post_adj)

    push!(post_networks, (
            fw_ID=string(fwid),
            Model=string(model_name),
            AdjacencyMatrix_burn = burn_adj_out,
            AdjacencyMatrix = post_adj_out,
            MetabolicClasses=out.met_class_post
        ); promote=true)
end


# --- 5. Save Simulation Summary ---
CSV.write("outputs/stability_metrics.csv", stability_metrics)
JLD2.save_object("networks/networks_END.jld2", post_networks)

