using LinearAlgebra # For matrix operations (e.g., transpose)
using Plots # For plotting the result

# --- 1. Parameters and Initial Condition ---
T = 50 # Years to simulate
x0 = [1.0; 0.0] # Initial distribution: [Share_Low; Share_High]

# --- 2. Simulation Function ---
function simulate_income_share(rho, T, x0)
    # Define the Transition Matrix P (from row to column)
    P = [ rho 1-rho; 1-rho rho ]
    
    # Store the shares (Low Income is the first element, index 1)
    share_L = zeros(T + 1) 
    xt = x0 
    share_L[1] = xt[1] 

    # Loop to simulate the distribution's evolution
    for t = 1:T
        # Calculate next period's distribution: x_t+1 = P' * x_t
        xt = transpose(P) * xt 
        share_L[t+1] = xt[1]
    end
    return share_L
end

# --- 1(a) Simulation for ρ = 0.8 ---
rho_low = 0.8
share_L_low_p = simulate_income_share(rho_low, T, x0)

# --- 1(b) Simulation for ρ = 0.99 ---
rho_high = 0.99
share_L_high_p = simulate_income_share(rho_high, T, x0)

# --- 3. Plotting the Results ---
time_vector = 0:T

plot(time_vector, share_L_low_p, 
    label="share low income (ρ=0.8)", 
    lw=2, 
    xlabel="year", 
    ylabel="Share in Low Income State",
    title="Convergence to Stationary Distribution",
    legend=:topright
)

plot!(time_vector, share_L_high_p, 
    label="share low income with high persistence (ρ=0.99)", 
    lw=2
)

hline!([0.5], line=:dash, color=:black, label="Stationary Share (0.5)")
# display(current())