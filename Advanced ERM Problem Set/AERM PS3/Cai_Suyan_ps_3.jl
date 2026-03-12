
# Exercise 1: Income Risk

## Q2: Stationary Distribution 
using LinearAlgebra

function rouwenhorst(rho, N, sigma_sq_u)
    # 1. Calculate parameters and grid width
    p = (1 + rho) / 2
    sigma_y = sqrt(sigma_sq_u / (1 - rho^2)) # Stationary variance formula 
    psi = sqrt(N - 1) * sigma_y
    
    # 2. Base Case N=2
    Pi = [p 1-p; 1-p p]
    
    # 3. Recursive Construction [cite: 291]
    for n in 3:N
        Pi_old = Pi
        Pi = zeros(n, n)
        
        Pi[1:n-1, 1:n-1] .+= p * Pi_old
        Pi[1:n-1, 2:n]   .+= (1-p) * Pi_old
        Pi[2:n, 1:n-1]   .+= (1-p) * Pi_old
        Pi[2:n, 2:n]     .+= p * Pi_old
        
        # Normalize middle rows [cite: 293]
        Pi[2:n-1, :] ./= 2
    end
    
    # 4. Create the log income grid (y)
    grid = collect(range(-psi, psi, length=N))
    
    return Pi, grid
end

function compute_stationary_distribution(Pi)
    # The stationary distribution is the left eigenvector corresponding to eigenvalue 1
    # of the transposed matrix Pi'
    
    # Eigenvalue decomposition of Pi'
    F = eigen(Pi') 
    
    # Find the eigenvector corresponding to eigenvalue 1 (or the one closest to 1)
    # We use a tolerance for comparison
    eig_val_1_index = argmin(abs.(F.values .- 1.0))
    
    # The eigenvector is complex-valued, take the real part
    pi_raw = real.(F.vectors[:, eig_val_1_index])
    
    # Normalize the eigenvector to sum to 1
    pi_stationary = pi_raw ./ sum(pi_raw)
    
    return pi_stationary
end

# RETURNING TO CALCULATION ---
rho = 0.9
sigma_sq_u = 0.015
Ny = 7

Pi_matrix, y_grid = rouwenhorst(rho, Ny, sigma_sq_u)
pi_dist = compute_stationary_distribution(Pi_matrix)

# Print results 
println("Log Income Grid (y):")
println(y_grid)
println("\nStationary Distribution (π):")
println(pi_dist)

# PLOTTING ---

using Plots
using Statistics # Needed for exp() function

# Given data from the previous calculation
# y_grid = [-0.6882, -0.4588, -0.2294, 0.0, 0.2294, 0.4588, 0.6882]
# pi_dist = [0.0156, 0.0938, 0.2344, 0.3125, 0.2344, 0.0938, 0.0156]

# 1. Calculate the Income Level (Y = e^y)
income_level = exp.(y_grid)

# 2. Generate the plot
# You can use a bar plot or a scatter plot with vertical lines (sticks)
# A bar plot is usually clearer for discrete distributions.
p = bar(income_level, pi_dist, 
    label="Stationary Distribution",
    title="Stationary Distribution over Income Level (eʸ)",
    xlabel="Income Level (Y = eʸ)",
    ylabel="Stationary Probability (π)",
    legend=false,
    bar_width=0.15, # Adjust width as needed for visual clarity
    xticks=income_level, # Ensure ticks are exactly at the income levels
    xrotation=45, # Rotate x-axis labels for readability
    bottom_margin=80, # increase bottom margin so xlabel isn't cut off (pixels)
    size=(900,500) # explicit size helps avoid cropping on save
)

# Save the figure to the local `figures/` folder (create folder if needed)
try
    mkpath("figures")
    savefig(p, "figures/stationary_income.png")
    println("Saved plot to figures/stationary_income.png")
catch e
    @warn "Could not save figure" exception=(e,catch_backtrace())
end






## Q3: Simulate Income Process
using LinearAlgebra
using Random
using Statistics

# --- 1. ROUWENHORST FUNCTION (from Problem 1) ---
function rouwenhorst(rho, N, sigma_sq_u)
    p = (1 + rho) / 2
    sigma_y = sqrt(sigma_sq_u / (1 - rho^2))
    psi = sqrt(N - 1) * sigma_y
    
    Pi = [p 1-p; 1-p p]
    
    for n in 3:N
        Pi_old = Pi
        Pi = zeros(n, n)
        
        Pi[1:n-1, 1:n-1] .+= p * Pi_old
        Pi[1:n-1, 2:n]   .+= (1-p) * Pi_old
        Pi[2:n, 1:n-1]   .+= (1-p) * Pi_old
        Pi[2:n, 2:n]     .+= p * Pi_old
        
        Pi[2:n-1, :] ./= 2
    end
    
    grid = collect(range(-psi, psi, length=N))
    return Pi, grid
end

# --- 2. STATIONARY DISTRIBUTION (from Problem 2) ---
function compute_stationary_distribution(Pi)
    # Finds the left eigenvector corresponding to eigenvalue 1 (Pi' * pi = 1 * pi)
    F = eigen(Pi') 
    eig_val_1_index = argmin(abs.(F.values .- 1.0))
    pi_raw = real.(F.vectors[:, eig_val_1_index])
    pi_stationary = pi_raw ./ sum(pi_raw)
    return pi_stationary
end

# --- 3. MARKOV SIMULATION LOGIC (Problem 3) ---

function markov_simulation(Pi, y_grid, pi_dist, N_sim, T_sim, seed)
    
    # Set the seed for reproducibility
    Random.seed!(seed)
    
    N_y = length(y_grid)
    
    # A. Generate Random Uniform Draws
    # Note: Problem 3 requires T=10, which means 1 initial draw + 9 transition draws.
    # We use a T_sim = 10 matrix and treat the first column as the initial draw.
    rand_mat_N_T = rand(N_sim, T_sim) # The requested rand_mat_N_T (2000 x 10)
    
    # B. Compute Conditional CDF (CCDF) Matrix F from Pi
    # F[i, j] is P(next state <= j | current state = i)
    F = cumsum(Pi, dims=2)
    
    # C. Compute Unconditional CDF (UCDF) for initial draw from stationary dist
    UCDF = cumsum(pi_dist)
    
    # Initialize matrices to store state indices (1 to Ny) and log income values (y)
    state_mat = zeros(Int, N_sim, T_sim)
    income_mat = zeros(N_sim, T_sim)
    
    # --- STEP 1: INITIAL DRAW (t=1) from Stationary Distribution ---
    # The first column of rand_mat_N_T is used for the initial draw
    for i in 1:N_sim
        # Find the state index (j) where the UCDF first exceeds the random draw
        j_initial = findfirst(x -> rand_mat_N_T[i, 1] < x, UCDF)
        
        state_mat[i, 1] = j_initial
        income_mat[i, 1] = y_grid[j_initial]
    end
    
    # --- STEP 2: RECURSIVE TRANSITION (t=2 to T_sim) ---
    for t in 2:T_sim
        for i in 1:N_sim
            # Current state index
            j_current = state_mat[i, t-1]
            
            # Find the state index (j_next) where the CCDF for the current state (j_current)
            # first exceeds the random draw for the next period.
            j_next = findfirst(x -> rand_mat_N_T[i, t] < x, F[j_current, :])
            
            state_mat[i, t] = j_next
            income_mat[i, t] = y_grid[j_next]
        end
    end
    
    return income_mat
end


# --- EXECUTION ---
rho = 0.9
sigma_sq_u = 0.015
N_y = 7
N_sim = 2000
T_sim = 10
seed = 1776

# Compute inputs needed for simulation
Pi_matrix, y_grid = rouwenhorst(rho, N_y, sigma_sq_u)
pi_dist = compute_stationary_distribution(Pi_matrix)

# Run the simulation
income_simulation_matrix = markov_simulation(Pi_matrix, y_grid, pi_dist, N_sim, T_sim, seed)

# Print a portion of the results for verification
println("\n--- SIMULATION RESULTS (First 5 households, all 10 years) ---")
for i in 1:5
    println("Household $i: ", income_simulation_matrix[i, :])
end




## Q4: Pooled OLS of Income Persistence
function estimate_ols_persistence(income_mat)
    
    N_sim, T_sim = size(income_mat)
    
    # We use periods t=2 to T_sim as the dependent variable (Y) 
    # and periods t=1 to T_sim-1 as the independent variable (X)
    
    # 1. Reshape Data for Regression
    # Y = y_i,t (t=2 to 10)
    Y = vec(income_mat[:, 2:T_sim]) # Flattens the 2000x9 matrix into a 18000-element vector
    
    # X = y_i,t-1 (t=1 to 9)
    # We use a matrix for X to easily handle the X'X and X'Y calculations later.
    # Since we are running OLS without an intercept, the X matrix is just one column (the regressor).
    X = reshape(vec(income_mat[:, 1:T_sim-1]), :, 1) # 18000 x 1 matrix
    
    # 2. Pooled OLS Calculation (without intercept: y = beta*x)
    # The formula is beta_hat = (X'X)^-1 * X'Y
    
    # Calculate X'X (a 1x1 scalar in this case)
    XTX = X' * X
    
    # Calculate X'Y (a 1x1 scalar in this case)
    XTY = X' * Y
    
    # Calculate the OLS coefficient
    beta_hat = XTX \ XTY # Julia's backslash operator solves XTX * beta_hat = XTY
    
    return beta_hat[1]
end


# Example: Run the OLS (Requires the previous code's output to be in memory)
beta_ols = estimate_ols_persistence(income_simulation_matrix)
println("\n--- Problem 4: OLS Persistence Estimate ---")
println("Estimated Persistence (β̂): ", beta_ols)


## Q5: Standard Error and Confidence Interval

function calculate_ols_se_and_ci(income_mat, rho_true)
    
    N_sim, T_sim = size(income_mat)
    
    # Data reshape (as in Problem 4)
    Y = vec(income_mat[:, 2:T_sim])
    X = reshape(vec(income_mat[:, 1:T_sim-1]), :, 1)
    N_obs = length(Y) # 18000 observations
    
    # 1. OLS Coefficient (from Problem 4 result)
    XTX = X' * X
    XTY = X' * Y
    beta_hat = XTX \ XTY
    
    # 2. Residuals (u_hat)
    u_hat = Y - X * beta_hat
    
    # 3. Sum of Squared Residuals (u_hat' * u_hat)
    u_hat_sq_sum = u_hat' * u_hat
    
    # 4. Variance of the Estimator: Var(beta_hat) = (X'X)^-1 * (u'u / N_obs)
    # Note: The problem asks for variance using the formula Var(beta) = (X'X)^-1 * (u'u / N).
    # We must be careful with the dimensions of (X'X)^-1. Since X'X is a 1x1 scalar, 
    # (X'X)^-1 is simply 1 / (X'X).
    Var_beta = (1 / XTX[1]) * (u_hat_sq_sum[1] / N_obs)
    
    # 5. Standard Error (SE)
    SE_beta = sqrt(Var_beta)
    
    # 6. 95% Confidence Interval (using Z-score of 1.96 for large N)
    z_score = 1.96
    lower_ci = beta_hat[1] - z_score * SE_beta
    upper_ci = beta_hat[1] + z_score * SE_beta
    
    # 7. Check if true rho is in CI
    rho_in_ci = (lower_ci <= rho_true <= upper_ci)
    
    return beta_hat[1], SE_beta, lower_ci, upper_ci, rho_in_ci
end

# --- EXECUTION ---
rho_true = 0.9

# Run the calculation (Requires the previous code's output to be in memory)
beta_ols, SE_beta, lower_ci, upper_ci, rho_in_ci = calculate_ols_se_and_ci(income_simulation_matrix, rho_true)

println("\n--- Problem 5: Standard Error and CI ---")
println("Estimated SE(β̂): ", SE_beta)
println("95% CI: [", lower_ci, ", ", upper_ci, "]")
println("True ρ=0.9 is in CI: ", rho_in_ci)

# Print results to be inserted into LaTeX
println("\n--- LaTeX INSERT VALUES ---")
println("SE(β̂): ", round(SE_beta, digits=5))
println("Lower CI: ", round(lower_ci, digits=4))
println("Upper CI: ", round(upper_ci, digits=4))



# Exercise 2: Income Shocks and Consumption Response
## Q1: Solve an Incomplete Markets Model 
using LinearAlgebra
using Interpolations # For linear interpolation
using Statistics

# --- 1. ROUWENHORST FUNCTION (Adapted for N_y=5) ---
function rouwenhorst(rho, N, sigma_sq_u)
    p = (1 + rho) / 2
    sigma_y = sqrt(sigma_sq_u / (1 - rho^2))
    psi = sqrt(N - 1) * sigma_y
    Pi = [p 1-p; 1-p p]
    
    for n in 3:N
        Pi_old = Pi
        Pi = zeros(n, n)
        Pi[1:n-1, 1:n-1] .+= p * Pi_old
        Pi[1:n-1, 2:n]   .+= (1-p) * Pi_old
        Pi[2:n, 1:n-1]   .+= (1-p) * Pi_old
        Pi[2:n, 2:n]     .+= p * Pi_old
        Pi[2:n-1, :] ./= 2
    end
    grid = collect(range(-psi, psi, length=N))
    return Pi, grid
end

# --- 2. MODEL PARAMETERS ---
const RHO = 0.98     # Income persistence
const SIGMA_U_SQ = 0.015 # Shock variance
const NY = 5         # Income states

const GAMMA = 2.0    # CRRA coefficient
const BETA = 0.95    # Discount factor
const R = 0.03       # Interest rate
const W = 1.0        # Wage
const R_GROSS = 1.0 + R

const NB = 100       # Number of asset grid points
const B_LOWER = 0.0  # Borrowing constraint (B >= 0)

# --- 3. UTILITY AND EXPECTED VALUE FUNCTION ---
function u(c)
    # CRRA utility: u(c) = c^(1-gamma) / (1-gamma)
    return (c^(1 - GAMMA) - 1) / (1 - GAMMA)
end

function expected_value(V_next, Pi_row)
    # E[V(B', y')] = Pi_row * V_next, where V_next is a column vector
    return Pi_row' * V_next
end

# --- 4. ASSET GRID SETUP ---
function setup_grids(B_upper_multiplier=15.0)
    # 4.1 Income Grid
    Pi, y_grid = rouwenhorst(RHO, NY, SIGMA_U_SQ)
    Y_grid = exp.(y_grid) # Income in levels

    # 4.2 Asset Grid
    Y_max = maximum(Y_grid)
    B_UPPER = B_upper_multiplier * Y_max # Following the rule: 15 * max_income (approx 51.3)
    
    # Log-linear spacing as requested: exp(LinRange(0, Log(B_upper+1), N)) - 1
    B_grid_log = collect(LinRange(0, log(B_UPPER + 1), NB))
    B_grid = exp.(B_grid_log) .- 1
    
    # 4.3 Policy Interpolation Grid (used for interpolation, B_PRIME_grid)
    B_prime_grid = copy(B_grid)
    
    return Pi, Y_grid, B_grid, B_prime_grid
end

# --- 5. VALUE FUNCTION ITERATION (VFI) ---
function solve_model_vfi(Pi, Y_grid, B_grid, B_prime_grid; max_iter=1000, tol=1e-6)
    
    NB_prime = length(B_prime_grid)
    NB_grid = length(B_grid)
    NY_grid = length(Y_grid)
    
    # Initialize Value Function V and Policy Function B_prime (B')
    V = zeros(NB_grid, NY_grid)
    V_new = zeros(NB_grid, NY_grid)
    B_prime = zeros(NB_grid, NY_grid) # Optimal next-period asset choice
    
    # Pre-calculate the Continuation Value (EV) term for all possible B_prime choices
    # This stores the expected utility from V(B', y') for all B' on the grid
    EV_matrix = zeros(NB_prime, NY_grid)
    
    # VFI Loop
    for iter in 1:max_iter
        
        # 1. Update Expected Value (EV) Matrix using current V
        for y_next_state in 1:NY_grid
            # Interpolate V for the next period over the B' grid
            # V[:, y_next_state] is the value function for a specific y' state across all B' assets.
            V_interpolator = LinearInterpolation(B_grid, V[:, y_next_state], extrapolation_bc=Line())
            EV_matrix[:, y_next_state] = V_interpolator.(B_prime_grid)
        end
        
        # 2. Calculate Expected Continuation Value (E_V) for each Pi row
        # E_V_at_B_prime_y is E[V(B', y') | y]
        # Size: NB_prime x NY_grid
        E_V_at_B_prime_y = EV_matrix * Pi' # Matrix multiplication: (NB' x NY) * (NY x NY) -> (NB' x NY)
        
        # 3. Bellman Equation Maximization
        V_max_diff = 0.0
        
        for y_state in 1:NY_grid # Loop over current income state (y)
            for b_idx in 1:NB_grid # Loop over current asset grid point (B)
                
                B = B_grid[b_idx]
                Y = Y_grid[y_state]
                
                # Current resources (Cash on Hand, CoH)
                CoH = R_GROSS * B + W * Y 
                
                max_val = -Inf
                best_B_prime = 0.0
                
                # Loop over possible next-period asset choices (B')
                for b_prime_idx in 1:NB_prime
                    B_prime_choice = B_prime_grid[b_prime_idx]
                    
                    # Ensure B' <= CoH (no negative consumption)
                    if B_prime_choice <= CoH 
                        c = CoH - B_prime_choice
                        
                        # Current Utility + Discounted Expected Future Value
                        val = u(c) + BETA * E_V_at_B_prime_y[b_prime_idx, y_state]
                        
                        if val > max_val
                            max_val = val
                            best_B_prime = B_prime_choice
                        end
                    else
                        # Since B_prime_grid is increasing, we can break early
                        break 
                    end
                end # end loop over B_prime_idx

                V_new[b_idx, y_state] = max_val
                B_prime[b_idx, y_state] = best_B_prime
            end # end loop over b_idx
        end # end loop over y_state
        
        # 4. Check Convergence
        V_max_diff = maximum(abs.(V_new .- V))
        if V_max_diff < tol
            println("VFI converged in $iter iterations.")
            break
        end
        
        V = V_new # Update V
        
        if iter == max_iter
            println("VFI failed to converge after $max_iter iterations (Max Diff: $V_max_diff).")
        end
    end # end VFI loop
    
    # Policy function c(B, y) is derived from B'(B, y)
    C_policy = R_GROSS * B_grid' .+ W * Y_grid' .- B_prime
    
    return B_prime, C_policy
end

# --- EXECUTION ---
Pi, Y_grid, B_grid, B_prime_policy = setup_grids()

# Solve the model (this step may take a significant amount of time to run)
# B_policy is the optimal B'(B, y), C_policy is the optimal c(B, y)
B_policy, C_policy = solve_model_vfi(Pi, Y_grid, B_grid, B_prime_policy)

# Print a summary of the resulting policy function (e.g., for the lowest income state)
println("\n--- SOLVED POLICY FUNCTION SUMMARY (Lowest Income State: Y=$(round(Y_grid[1], digits=3))) ---")
println("B_grid (first 5): ", round.(B_grid[1:5], digits=4))
println("C_policy (first 5): ", round.(C_policy[1:5, 1], digits=4))
println("B_prime_policy (first 5): ", round.(B_policy[1:5, 1], digits=4))