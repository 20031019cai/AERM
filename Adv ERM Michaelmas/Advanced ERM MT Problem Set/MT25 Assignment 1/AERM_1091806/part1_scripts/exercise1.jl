# Exam MT 2025 Assignment 1
# Exercise 1 Part 1

# === Set-up ===

# This is for saving figures to the results folder
const RESULTS_DIR = "/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results"

# Load packages
using Pkg
using Plots
using Random
using Statistics
using LinearAlgebra

# ==============================================================================
# Question 1
# ==============================================================================

# === Set-up ===
# Income states

    yL = 1.0
    
    yH = 2.0

# Values of persistence to consider

    ρ_vals = [0.0, 0.05, 0.2, 0.5, 0.7, 0.95, 1.0]

# Time horizon
    T = 100

# Initial distribution: P(y = yL) = 0.2
    π0 = [0.2, 0.8]  # [P(yL), P(yH)]using Pkg

# Transition matrix for symmetric persistence ρ:
#       yL      yH
# yL [  ρ    1-ρ  ]
# yH [ 1-ρ    ρ   ]
    function transition_matrix(ρ)
        return [ρ (1-ρ); (1-ρ) ρ]
    end

# Evolve distribution forward: πₜ₊₁ = πₜ' * P
# (row vector times transition matrix)
    function evolve_distribution(π0, P, T)
        shares_low = zeros(T + 1)
        shares_low[1] = π0[1]
        πt = copy(π0)
        for t in 1:T
            πt = πt' * P |> vec
            shares_low[t + 1] = πt[1]
        end
        return shares_low
    end

# Compute and plot
    Pkg.add("Plots")
    plt = plot(
        xlabel = "Year",
        ylabel = "Share of low-income households",
        title = "Evolution of low-income share by persistence ρ",
        legend = :right,
        size = (800, 500),
        linewidth = 2
    )

    for ρ in ρ_vals
        P = transition_matrix(ρ)
        shares = evolve_distribution(π0, P, T)
        plot!(plt, 0:T, shares, label = "ρ = $ρ")
    end

# Add reference line at 0.5

    hline!(plt, [0.5], color = :black, linestyle = :dash, label = "0.5", alpha = 0.5)

# ******* Save figure  ****** [TO DELETE BEFORE SUBMISSION]
    mkpath("/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results/")
    savefig(plt, "/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results/fig_ex1_q1.pdf")
    display(plt)

# ==============================================================================
# Question 2
# ==============================================================================

# At the stationary distribution (50/50), Pr(switch) = 1-ρ for each type
# Δlog(y) ∈ {0, log(2), -log(2)}

    Δ = log(2)

    println("ρ \t\t Var(Δlog y)")

    println("—"^30)

    for ρ in ρ_vals
        # Stationary shares
            p_L = (ρ == 1.0) ? 0.5 : 0.5  # symmetric chain

            p_H = 1 - p_L
        
        # Joint probabilities of transitions

            p_LL = p_L * ρ

            p_LH = p_L * (1 - ρ)
            
            p_HL = p_H * (1 - ρ)
            
            p_HH = p_H * ρ
            
        # E[Δlog y] and E[(Δlog y)²]

            mean_Δ = p_LH * Δ + p_HL * (-Δ)

            mean_Δ_sq = p_LH * Δ^2 + p_HL * Δ^2
            
            var_Δ = mean_Δ_sq - mean_Δ^2

        println("ρ = $ρ \t\t $(round(var_Δ, digits=5))")
    end


# ==============================================================================
# Question 3
# ==============================================================================

# Set seed

    Random.seed!(1776)

# Generate random numbers for 1000 households and 8 years

    rand_mat = rand(1000, 8)

# Set parameters

    N = 1000
    T_sim = 8
    ρ_sim = 0.95
    

# Transition matrix for ρ = 0.95

    P_sim = transition_matrix(ρ_sim)

# Allocate state matrix: 1 = low income, 2 = high income
# Use column 1 of rand_mat to assign initial condition:
# ~20% start in low state, ~80% start in high state
    states = fill(2, N, T_sim)
    
    states[rand_mat[:, 1] .< 0.2, 1] .= 1

# Simulate forward using columns 2:T_sim of rand_mat
# For each household, compare rand draw to cumulative probability of staying in current state
# If rand < P[s,1] → go to state 1 (low), otherwise → state 2 (high)

    for t in 2:T_sim
        for i in 1:N

            # Get current state
            s = states[i, t-1]

            # Update state based on random draw and transition probabilities
            states[i, t] = rand_mat[i, t] < P_sim[s, 1] ? 1 : 2
        end
    end

# Compute share of high types in simulation at each period
   
    share_high_sim = vec(mean(states .== 2, dims=1))

# Compute share of high types from the continuum (theoretical) exercise
# Reuse the evolve_distribution function from Q1 with ρ = 0.95
    
    P95 = transition_matrix(0.95)
    π_cont = [0.2, 0.8]  # initial distribution [P(low), P(high)]
    shares_low_cont = evolve_distribution(π_cont, P95, T_sim - 1)
    share_high_cont = 1.0 .- shares_low_cont  # convert to high-type share

# Plot comparison

    plt3 = plot(1:T_sim, share_high_sim,

        marker=:circle, linewidth=2, label="1000 households",
        xlabel="Year", ylabel="Share of high-income households",
        title="Simulated vs. Continuum: Share of High Types (ρ=0.95)",
        size=(800, 500))

    plot!(plt3, 1:T_sim, share_high_cont,
        marker=:circle, linestyle=:dash, linewidth=2, label="Continuum")

# ******* Save figure  ****** [TO DELETE BEFORE SUBMISSION]

    mkpath("/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results/")
    savefig(plt3, "/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results/fig_ex1_q3.pdf")
    display(plt3)

# ==============================================================================
# Question 4
# ==============================================================================

# Convert states to income levels (1 = yL = 1, 2 = yH = 2)

    Y_panel = float.(states)
    Y_panel[states .== 1] .= yL
    Y_panel[states .== 2] .= yH

# Dependent variable y_{i,t} for t=2,...,8 and lag y_{i,t-1} for t=1,...,7
    
    Y_curr = Y_panel[:, 2:end]
    Y_lag  = Y_panel[:, 1:end-1]

# Stack into vectors

    y_vec     = vec(Y_curr)
    y_lag_vec = vec(Y_lag)

# --- Model 1: Pooled OLS ---

    X_ols = hcat(ones(length(y_lag_vec)), y_lag_vec)
    β_ols = X_ols \ y_vec
    
    println("OLS:  α = $(round(β_ols[1], digits=4)),  β = $(round(β_ols[2], digits=4))")

# --- Model 2: Fixed Effects (within transformation) ---

    y_curr_mean = mean(Y_curr, dims=2)
    y_lag_mean  = mean(Y_lag, dims=2)

    Y_curr_dm = Y_curr .- y_curr_mean
    Y_lag_dm  = Y_lag  .- y_lag_mean

    y_dm     = vec(Y_curr_dm)
    y_lag_dm = vec(Y_lag_dm)

    β_fe = dot(y_lag_dm, y_dm) / dot(y_lag_dm, y_lag_dm)

    println("FE:   β = $(round(β_fe, digits=4))")
    println("True β = 2ρ - 1 = $(2*0.95 - 1)")
    println("OLS: α = $(round(β_ols[1], digits=4)), β = $(round(β_ols[2], digits=4))
")












