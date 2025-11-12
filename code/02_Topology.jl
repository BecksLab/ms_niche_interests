#=
-------------------------------------------------
02_Topology.jl
Calculating structure/topology for generated networks.
-------------------------------------------------

=#
# Load the Pkg package manager.
using Pkg
# Activate the current project environment to ensure consistent package versions.
Pkg.activate(".")

# --- 1. Load Dependencies ---
using CSV
using DataFrames
using JLD2
using LinearAlgebra
using Plots
using SpeciesInteractionNetworks

# --- 2. Load All Code ---
include("lib/internals.jl");

# --- 3. Import networks .jld2 object ---

networks = load_object("data/outputs/network_test_verified_seed_42_29-10-2025.jld2")

# --- 4. Convert networks to `SpeciesInteractionNetworks` networks ---

# we will also append this to the network object as a col.
networks.InteractionNetwork = build_network.(networks.AdjacencyMatrix)

# --- 5. Get topology ---

# create df to store outputs

topology = DataFrame(
    model = String[],
    richness = Int64[],
    connectance = Float64[],
    complexity = Float64[],
    trophic_level = Float64[],
    generality = Float64[],
    vulnerability = Float64[],
    S1 = Float64[],
    S2 = Float64[],
    S4 = Float64[],
    S5 = Float64[],
);

# note for simplicity we use a wrapper function

for i in 1:nrow(networks)

    d = network_summary(networks.InteractionNetwork[i])

    d[:model] = networks.Model[i]

    push!(topology, d)

end

# write summaries as .csv
CSV.write("data/outputs/topology.csv", topology)

# --- 6. Linear Discriminant Analysis (LDA) ---

# prep data
x = Matrix(topology[:, 3:end])
x_labels = convert(Vector{String}, topology[:, :model])

# Compute overall mean
μ = mean(x, dims=1)

classes = unique(x_labels)
n_features = size(x, 2)

# Within-class scatter (SW) and between-class scatter (SB)
SW = zeros(n_features, n_features)
SB = zeros(n_features, n_features)

for c in classes
    Xc = x[x_labels .== c, :]
    μc = mean(Xc, dims=1)
    SW += cov(Xc) * (size(Xc, 1) - 1)
    n_c = size(Xc, 1)
    mean_diff = (μc - μ)'
    SB += n_c * (mean_diff * mean_diff')
end

# Solve the generalized eigenvalue problem: inv(SW)*SB
eigenvals, eigenvecs = eigen(Symmetric(inv(SW) * SB))

# Sort eigenvectors by eigenvalues (descending)
sorted_idx = sortperm(eigenvals, rev=true)
W = eigenvecs[:, sorted_idx[1:2]]   # top 2 discriminants

# Project data
X_lda = x * W

# Plot

scatter()
for c in classes
    idx = findall(x_labels .== c)
    scatter!(X_lda[idx, 1], X_lda[idx, 2], label=c)
end
xlabel!("LD1")
ylabel!("LD2")
title!("LDA Projection")