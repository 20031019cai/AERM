# Exam MT 2025 Assignment 1
# Part 1 Exercise 2

# === Set up ===

    # Load packages

        using LinearAlgebra
        using Statistics
        using Plots
        using Printf

    # Set up paths

        const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

    # Load code files for functions

        include("ex2_parameters.jl")
        include("ex2_vfi.jl")
        include("ex2_kfe.jl")

    # Set up parameters and ρ values    

        p = Params()
        ρ_vals = [0.0, 0.05, 0.2, 0.5, 0.7, 0.95, 1.0]

# ============================================================
# Question 1
# ============================================================

asset_supply = Dict{Float64, Float64}()
solutions = Dict{Float64, Any}()

for ρ in ρ_vals
    println("\n=== Solving for ρ = $ρ ===")
    sol = solve_household(p, ρ)
    solutions[ρ] = sol

    if ρ == 1.0
        # ρ=1: absorbing states, no unique stationary dist
        # Under 50/50 split, both types hit borrowing constraint → B = -0.5
        asset_supply[ρ] = -0.5
        println("ρ=1.0: absorbing states → B = -0.5 (by argument)")
    else
        μ = stationary_dist(sol, ρ, p)
        B = agg_asset_supply(μ, sol)
        asset_supply[ρ] = B
        println("Asset supply: B = $B")
    end
end

# Print table
println("\n===== Q1: Asset Supply by ρ =====")
println("ρ \t\t Asset Supply")
println("—"^35)
for ρ in ρ_vals
    @printf("%.2f \t\t %.6f\n", ρ, asset_supply[ρ])
end



# ============================================================
# Question 2a: Counterfactual asset supply
# Keep c = c^{0.5} but use income process for ρ=0 and ρ=0.95
# ============================================================

# Benchmark solution: ρ = 0.5
sol_benchmark = solutions[0.5]

# Counterfactual: use benchmark policy but different ρ for income transitions
B_cf = Dict{Float64, Float64}()
for ρ_cf in [0.0, 0.95]
    # Use the benchmark solution's policy_idx, but compute stationary dist under ρ_cf
    μ_cf = stationary_dist(sol_benchmark, ρ_cf, p)
    B_cf[ρ_cf] = agg_asset_supply(μ_cf, sol_benchmark)
    println("B(ρ=$ρ_cf, c=c^0.5) = $(B_cf[ρ_cf])")
end



# ============================================================
# Question 2b: Decomposition table
# ============================================================

B_benchmark = asset_supply[0.5]

# ρ = 0 decomposition
T_total_0 = asset_supply[0.0] - B_benchmark
T1_0 = asset_supply[0.0] - B_cf[0.0]
T2_0 = B_cf[0.0] - B_benchmark

# ρ = 0.95 decomposition
T_total_95 = asset_supply[0.95] - B_benchmark
T1_95 = asset_supply[0.95] - B_cf[0.95]
T2_95 = B_cf[0.95] - B_benchmark

# Print decomposition table
println("\n===== Q2b: Decomposition =====")
println("ρ \t T_total \t T1 (%) \t\t T2 (%)")
println("—"^55)
@printf("0.00 \t %.4f \t\t %.1f \t\t %.1f\n", T_total_0, 100*T1_0/T_total_0, 100*T2_0/T_total_0)
@printf("0.95 \t %.4f \t\t %.1f \t\t %.1f\n", T_total_95, 100*T1_95/T_total_95, 100*T2_95/T_total_95)

