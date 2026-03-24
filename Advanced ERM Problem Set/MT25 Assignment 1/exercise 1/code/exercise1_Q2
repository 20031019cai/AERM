# Advanced ERM Problem Set, MT25 Examined Assignment 1

# Exercise 1

# Q2: Variance of change in log labour income

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