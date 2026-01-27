
# Solving systems of equations numerically 
import Ipopt
using JuMP 

# Doing statistics and models
using Distances
using Distributions 
using LinearAlgebra 
using ProgressMeter
using Random
using SparseArrays
using Statistics
using StatsBase
using Turing

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