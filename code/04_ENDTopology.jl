#=
-------------------------------------------------
04_ENDTopology.jl
Calculating equilibrium structure/topology.
-------------------------------------------------
=#

using Pkg
Pkg.activate(".")

Pkg.add(url="https://github.com/PoisotLab/SpeciesInteractionNetworks.jl.git")

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using LinearAlgebra
using Plots
using SpeciesInteractionNetworks
using Graphs 

# --- 2. Load All Code ---

include(joinpath("lib", "internals.jl"));

# --- 3. Import networks .jld2 object ---
date_str = "11-06-2026"
path = joinpath(BASE_DIR, "data", "outputs", "END_final_adj_$(date_str).jld2")
networks = load_object(path)

# --- 4. Convert networks to `SpeciesInteractionNetworks` ---
networks.InteractionNetwork = [
    ismissing(A) ? missing : build_network(Int.(A))
    for A in networks.alive_connected_A
]

# --- 5. Custom Equilibrium Topology Function ---
function calculate_equ_custom_metrics(A::Matrix{Int})
    S = size(A, 1)
    L = sum(A)
    conn = S > 0 ? L / (S^2) : 0.0
    
    in_deg = vec(sum(A, dims=1))
    out_deg = vec(sum(A, dims=2))
    
    isolated = sum((in_deg .== 0) .& (out_deg .== 0))
    prop_basal = S > 0 ? sum(in_deg .== 0) / S : 0.0
    
    g = SimpleDiGraph(A)
    cycles = simplecycles(g)
    closed_loops = length(cycles)
    
    # Illogical tracking at equilibrium requires trait mapping logic, 
    # omitted here as dynamic indices shift post-extinction.
    return (S, conn, prop_basal, isolated, missing, closed_loops)
end

# --- 6. Get topology ---
topology = DataFrame(
    fw_ID = String[],
    model = String[],
    S_equ = Int[],
    connectance_equ = Float64[],
    prop_basal_equ = Float64[],
    mx_tl_equ = Float64[],
    mean_FCL_equ = Float64[],
    closed_loops_equ = Int[],
    isolated_equ = Int[],
    illogical_equ = Union{Int, Missing}[],
);


for i in 1:nrow(networks)

    N = networks.InteractionNetwork[i]
    A_mat = networks.alive_connected_A[i]

    if ismissing(N) || ismissing(A_mat)
        continue
    end

    d_tl = NaN
    d_fcl = NaN
    try
        d = network_summary(N)
        d_tl = d[:trophic_level]
        d_fcl = d[:ChLen]
    catch
        @warn "network_summary failed for $(networks.fw_ID[i]). Assigning NaN."
    end
        
    S_eq, conn_eq, p_basal_eq, iso_eq, illog_eq, loops_eq = calculate_equ_custom_metrics(Int.(A_mat))

    push!(topology, (
        fw_ID = networks.fw_ID[i],
        model = networks.model[i],
        S_equ = S_eq,
        connectance_equ = conn_eq,
        prop_basal_equ = p_basal_eq,
        mx_tl_equ = d_tl, 
        mean_FCL_equ = d_fcl,      
        closed_loops_equ = loops_eq,
        isolated_equ = iso_eq,
        illogical_equ = illog_eq        
    ))
end

outpath = joinpath("data", "outputs", "END_topology_equ_$(date_str).csv")
CSV.write(outpath, topology)