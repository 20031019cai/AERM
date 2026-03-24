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


# Question 1: Stationary income distribution dynamics

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

# Save figure
mkpath("/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results/")
savefig(plt, "/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results/fig_ex1_q1.pdf")
display(plt)



# ============================================================
# Q3: Simulating income process for 1000 households
# ============================================================
Random.seed!(1776)
rand_mat = rand(1000, 8)

N = 1000
T_sim = 8
ρ = 0.95

# Initial states: 20% low (1), 80% high (2)
# 1 = low, 2 = high
states = ones(Int, N, T_sim)
states[1:Int(0.2*N), 1] .= 1      # first 200 are low
states[Int(0.2*N)+1:end, 1] .= 2   # rest are high

P = transition_matrix(ρ)

# Simulate forward: if rand < P[current_state, 1], go to state 1, else state 2
for t in 2:T_sim
    for i in 1:N
        s = states[i, t-1]
        states[i, t] = rand_mat[i, t] < P[s, 1] ? 1 : 2
    end
end

# Share of high types in simulation
share_high_sim = [mean(states[:, t] .== 2) for t in 1:T_sim]

# Continuum evolution with ρ = 0.95
P95 = transition_matrix(0.95)
π_cont = [0.2, 0.8]  # [low, high]
share_high_cont = zeros(T_sim)
share_high_cont[1] = π_cont[2]
πt = copy(π_cont)
for t in 2:T_sim
    πt = vec(πt' * P95)
    share_high_cont[t] = πt[2]
end

# Plot comparison
plt3 = plot(1:T_sim, share_high_sim, label="Simulation (N=1000)", marker=:circle,
            xlabel="Year", ylabel="Share of high-income households",
            title="Simulated vs. Continuum: Share of High Types (ρ=0.95)",
            size=(800, 500), linewidth=2)
plot!(plt3, 1:T_sim, share_high_cont, label="Continuum", linestyle=:dash, linewidth=2)


mkpath("/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results/")
savefig(plt3, "/Users/caisuyan/Documents/GitHub/AERM/Advanced ERM Problem Set/MT25 Assignment 1/exercise 1/results/fig_ex1_q3.pdf")
display(plt3)


