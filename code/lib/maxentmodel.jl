
# Solving systems of equations numerically 
import Ipopt
using JuMP 

# Doing statistics and models
using LinearAlgebra 
using ProgressMeter
using Random
using StatsBase

"""
solve_lambdas(S::Int64, L::Int64) 
    S: number of species
    L: number of links
Returns the numerical approximation of the Lagrange multipliers for the joint degree distribution of maximum entropy
"""
function solve_lambdas(S::Int64, L::Int64) 
    model = Model(Ipopt.Optimizer) # create a non-linear model
    set_silent(model)
    
    @variable(model, x[1:2] >= 0)  # lambdas (x) are positive real numbers 
    @NLobjective(model, Max, x[1]) # could have been another objective function
    @NLobjective(model, Max, x[2]) # could have been another objective function
    @NLconstraint(model, sum(sum(k_in * exp(-x[1]*k_in - x[2]*k_out) for k_in = 0:S) for k_out = 0:S) / sum(sum(exp(-x[1]*k_in - x[2]*k_out) for k_in = 0:S) for k_out = 0:S) == L/S) # constraint for the average in-degree
    @NLconstraint(model, sum(sum(k_out * exp(-x[1]*k_in - x[2]*k_out) for k_in = 0:S) for k_out = 0:S) / sum(sum(exp(-x[1]*k_in - x[2]*k_out) for k_in = 0:S) for k_out = 0:S) == L/S) # constraint for the average out-degree
    optimize!(model) # find numerical solution
    return value.(x) 
end

"""
joint_degree_dist_maxent(S::Int64, L::Int64) 
    S: number of species
    L: number of links
Returns the joint degree distribution of maximum entropy given S and L
"""
function joint_degree_dist_maxent(S::Int64, L::Int64) 
    x = solve_lambdas(S, L) # approximate the Lagrange multipliers
    
    joint_degree_dist = zeros(S+1, S+1) # create matrix for the joint degree distribution (S+1 rows and columns because in and out-degrees go from 0 to S)

    Z = sum(sum(exp(-x[1]*k_in - x[2]*k_out) for k_in = 0:S) for k_out = 0:S) # compute partition function (denominator of the MaxEnt distribution)

    # compute probability for every combination of in and out-degrees
    for k_in in 0:S
        for k_out in 0:S
            num = exp(-x[1]*k_in - x[2]*k_out) 
            joint_degree_dist[k_in + 1, k_out + 1] = num/Z
        end
    end   
    return joint_degree_dist                                 
end

"""
simulate_degrees(JDD::Matrix{Float64})
    JDD: joint degree distribution 
Returns a simulated vector of in and out degrees using the joint degree distribution as weight
"""
function simulate_degrees(JDD::Matrix{Float64})
    S = size(JDD, 1) - 1 
    deg = findall(JDD .>= 0) 
    deg_samp = sample(deg, Weights(vec(JDD)), S, replace=true) 
    
    kin = [d[1] - 1 for d in deg_samp]  # Index 1 is the row (k_in)
    kout = [d[2] - 1 for d in deg_samp] # Index 2 is the col (k_out)
    
    # Validation: In a food web, sum(kin) must equal sum(kout)
    return (kin = kin, kout = kout)
end

using LinearAlgebra, Random, ProgressMeter

"""
svd_entropy(A::Matrix{Bool})
Calculates the Shannon entropy of the normalized singular values of the matrix.
This is the "complexity" metric the paper seeks to maximize.
"""
function svd_entropy(A::AbstractMatrix)
    σ = svdvals(Float64.(A))
    s_sum = sum(σ)
    s_sum == 0 && return 0.0
    
    ω = σ ./ s_sum
    # Calculate -sum(p log p), filtering out near-zero values
    return -sum(x * log(x) for x in ω if x > 1e-12)
end

"""
initial_matrix(kin::Vector{Int}, kout::Vector{Int})
Creates a starting binary matrix from degree sequences using a configuration-style approach.
"""
function initial_matrix(kin::Vector{Int}, kout::Vector{Int})
    S = length(kin)
    A = zeros(Bool, S, S)
    
    # Create "stubs" for links
    in_stubs = vcat([fill(i, kin[i]) for i in 1:S]...)
    out_stubs = vcat([fill(i, kout[i]) for i in 1:S]...)
    
    shuffle!(in_stubs)
    shuffle!(out_stubs)
    
    # Connect stubs (respecting the smaller sum to avoid errors)
    L = min(length(in_stubs), length(out_stubs))
    for i in 1:L
        u, v = out_stubs[i], in_stubs[i]
        # Avoid self-loops and multi-edges in the seed
        if u != v && !A[u, v]
            A[u, v] = true
        end
    end
    return A
end

"""
build_maxent_network(kin, kout; iterations=5000, T0=0.1)
The "Heuristic Type II" model from the paper. 
Optimizes the matrix structure while keeping the degree of every species fixed.
"""
function build_maxent_network(kin, kout; iterations=5000, T0=0.1)
    # Start with a random seed matching the degrees
    A = initial_matrix(kin, kout)
    curr_h = svd_entropy(A)
    
    T = T0
    alpha = 0.999 # Cooling rate
    best_A = copy(A)
    best_h = curr_h

    @showprogress 1 "Optimizing Network Structure..." for i in 1:iterations
        # 1. Propose a degree-preserving swap (Switching Move)
        edges = findall(A)
        length(edges) < 2 && break
        
        # Pick two random directed edges: (u -> v) and (x -> y)
        idx1, idx2 = rand(1:length(edges), 2)
        u, v = edges[idx1].I
        x, y = edges[idx2].I
        
        # 2. Check if swapping to (u -> y) and (x -> v) is valid
        # (Ensures no self-loops and no duplicate links)
        if u != y && x != v && !A[u, y] && !A[x, v]
            # Apply swap
            A[u, v] = false; A[x, y] = false
            A[u, y] = true; A[x, v] = true
            
            new_h = svd_entropy(A)
            Δh = new_h - curr_h
            
            # 3. Metropolis Criterion (Maximize Entropy)
            if Δh > 0 || rand() < exp(Δh / T)
                curr_h = new_h
                if curr_h > best_h
                    best_h = curr_h
                    best_A = copy(A)
                end
            else
                # Reject: Revert the swap
                A[u, v] = true; A[x, y] = true
                A[u, y] = false; A[x, v] = false
            end
        end
        
        # Cool down the temperature
        if i % 10 == 0; T *= alpha; end
    end
    
    return best_A
end