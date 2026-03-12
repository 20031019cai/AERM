# Exercise 1
### 2. Print
println("Hello, World.")


### 3. Define the function f(x) ###
using Plots
f(x) = x^2 * (1 - x)

# Define the range
x_range = -5.0:0.01:5.0 # Create a range of x values. Start: -5.0, step: 0.01 (increments), end: 5.0
y_values = f.(x_range)   # Apply the function to all x values

# Plot the function
plot(x_range, y_values,
     label = "f(x) = x²(1-x)",
     xlabel = "x",
     ylabel = "f(x)",
     title = "Plot of f(x)",
     legend = :topright)


### 4. Roots of the function f(x) ###
# The roots are at x = 0 and x = 1

### 5. Bisection Method Implementation ###
# Define the function (re-defined for clarity)
f(x) = x^2 * (1 - x)

function bisection(f, a, b, tol)
    # Check if a root is bracketed (f(a) and f(b) must have opposite signs)
    if f(a) * f(b) >= 0
        error("Bisection method fails. f(a) and f(b) must have opposite signs.")
    end

    c = (a + b) / 2 # Midpoint
    
    # Loop until the function value at the midpoint is below tolerance
    while abs(f(c)) > tol
        c = (a + b) / 2
        # Check sign change: if f(a) and f(c) have opposite signs, root is in [a, c]
        if f(a) * f(c) < 0
            b = c
        else
            # Otherwise, the root is in [c, b] (since f(b) must have the opposite sign of f(a))
            a = c
        end
    end
    return c
end

# Define a tolerance
tolerance = 1e-6

# --- Finding Root 2 (x=1) ---
# Bracketing interval: [0.5, 2.0]
# f(0.5) = 0.125 (Positive)
# f(2.0) = -4.0  (Negative)
# Bisection condition f(a)*f(b) < 0 is satisfied
root2_custom = bisection(f, 0.5, 2.0, tolerance)
println("Custom Bisection Root 2 (x=1): ", root2_custom)


### 6. Simulation and Seed Setting ###
# 6. (a) Explain what it means to set the seed
     # 1. What it is: setting the seed involves initializing the Pseudorandom Number Generator (PRNG) with a specific starting value

     # 2. Why it's useful: Setting the seed ensures reproducibility. If we run the code multiple times with the same seed, you will always get the exact same sequence of 20,000 "random" draws. This is essential for debugging, sharing results, and making your research verifiable

# 6. (b) Simulation code
using Random # Required for the `seed!` function

# Set the seed to 1776 for reproducibility
Random.seed!(1776)

# Simulate 20,000 draws from a standard uniform distribution (on [0, 1])
N_draws = 20000
uniform_draws = rand(N_draws)

# Display the first few draws to confirm:
println("First 5 uniform draws: ", uniform_draws[1:5])


### 7. Bernoulli Process Simulation

using Distributions # For Bernoulli distribution
using Plots       # For plotting histograms

# Define the function as suggested in the hint 
function simulate_bernoulli_means(N, p, R)
    # N: number of observations per simulation (40)
    # p: Bernoulli success probability (0.5)
    # R: number of repetitions (500)

    # 1. Define the Bernoulli distribution
    D = Bernoulli(p)

    # 2. Simulate R repetitions
    means = zeros(R) # Pre-allocate an array for R mean values

    for i in 1:R
        # Simulate N draws from the Bernoulli distribution (0 or 1)
        draws = rand(D, N)
        # Calculate the mean of the draws and store it
        means[i] = mean(draws)
    end
    return means
end

# Set parameters for the first run
N_observations = 40
p_success = 0.5
N_repetitions = 500

# Run the simulation
mean_distribution = simulate_bernoulli_means(N_observations, p_success, N_repetitions)

# Define the bins for the histogram 
# range [0, 1] with steps of 0.02
bin_edges = 0.0:0.02:1.0

# Plot the histogram
# We save the plot so we can easily update it for future exercises [cite: 34]
histogram(mean_distribution,
          bins = bin_edges,
          normalize = :probability, # Or :density, depending on preference
          legend = false,
          xlabel = "Mean of 40 Bernoulli Draws",
          ylabel = "Probability Density",
          title = "Histogram of Means (N=40, p=0.5, R=500)")

### 8. Repeat Simulation for varying p
using Random # We need this to reset the seed for each run
using Distributions
using Plots

# Re-define the parameters from step 7
N_observations = 40
N_repetitions = 500
P_values = [0.75, 0.9, 0.95, 0.99] # The new p values to test

# Redefine the function (or ensure it's loaded from step 7)
function simulate_bernoulli_means(N, p, R)
    D = Bernoulli(p)
    means = zeros(R)
    for i in 1:R
        # Use the Random package for simulation
        draws = rand(D, N)
        means[i] = sum(draws) / N
    end
    return means
end

# Define bins
bin_edges = 0.0:0.02:1.0  # from 0 to 1 with step 0.02

# Generate a combined plot (e.g., a 2x2 layout)
p_plots = []

for p in P_values
    # Note: The problem asks to use the "same uniform draws," which is tricky 
    # to enforce perfectly across runs, but setting the seed before *each* call 
    # ensures the sequence *starts* the same for each p, though the process differs.
    # The more standard interpretation for this kind of problem is to just
    # ensure overall reproducibility, which we do with the seed (step 6).

    # Simulate and calculate means
    mean_distribution = simulate_bernoulli_means(N_observations, p, N_repetitions)

    # Plot the histogram
    h = histogram(mean_distribution,
                  bins = bin_edges,
                  normalize = :probability,
                  legend = false,
                  title = "p = $(p)",
                  xlims = (0, 1))
    push!(p_plots, h)
end

# Combine the plots into a single figure
plot(p_plots..., layout = (2, 2), size = (800, 600),
     plot_title = "Distribution of Bernoulli Means (N=40, R=500)")

### 9. Large-Scale Simulation and Finer Grid

# Re-define parameters for the large experiment
N_observations_large = 20000     # Observations per simulation
p_success_large = 0.98  # New Bernoulli probability
N_repetitions_large = 100 # Number of repetitions

# Simulate and calculate means
@time mean_distribution_large = simulate_bernoulli_means(N_observations_large, p_success_large, N_repetitions_large)

# Define the finer bins for the histogram
# The means will be centered around p = 0.98. We need a range around this
# A reasonable range might be [0.95, 1.01] with a step of 0.002

bin_edges_fine = 0.95:0.002:1.01

# Plot the histogram, focusing on making it look nice
histogram(mean_distribution_large,
          bins = bin_edges_fine,
          normalize = :probability,
          legend = false,
          color = :blue,
          alpha = 0.7, # Transparency
          xlabel = "Mean of 20,000 Bernoulli Draws",
          ylabel = "Probability Density",
          title = "CLT Convergence: N=20000, p=0.98 (Fine Grid)",
          titlefontsize = 12,
          labelfontsize = 10,
          tickfontsize = 8,
          grid = true)

# Add a vertical line for the expected mean (p=0.98)
vline!([p_success_large], color = :red, linestyle = :dash, label = "Expected Mean")



# Exercise 2: Cake-eating problem 🎂
### 2. Solve the model computationally using value function iteration
# (a) Argue intuitively why this grid might lead to problems
     # The grid defined by LinRange(0, 10, 200) includes the state b=0
     # At b = 0, the only feasible consumption is c = 0 (since c =< 0 must hold). The standard utility function is {u(c) = log(c)
     # The problem is that consuming c = 0 yields u(0) - infty. This infinite negative value will contaminate the calculated value function, V(b), leading to numerical instability and overflow errors in the Value Function Iteration (VFI) process

# (b) Implement VFI with Bounded Utility and Plot Policy
using Plots        # For plotting the policy function
using Interpolations # For interpolating the value function V(b')

# Model Parameters
beta = 0.9           # Discount factor
c_min_bound = 0.01   # Consumption lower bound for utility fix

# Discretization Grid
N_b = 200
b_grid = LinRange(0.0, 10.0, N_b)

# --- 1. Define the Modified Utility Function ---
function u(c, c_min)
    # Utility is log(c), or log(c_min) if c is too small
    if c > c_min
        return log(c)
    else
        return log(c_min)
    end
end

# --- 2. VFI Initialization ---
tol = 1e-6           # Convergence tolerance
max_iter = 500       # Maximum iterations

V_old = zeros(N_b)
V_new = zeros(N_b)
policy = zeros(N_b) # Stores optimal consumption c(b)

# --- 3. VFI Main Loop ---
iter = 0
max_diff = 1.0 # Initialize max difference

while max_diff > tol && iter < max_iter
    
    # Create an interpolator for V_old
    V_old_interp = LinearInterpolation(b_grid, V_old, extrapolation_bc = Interpolations.Line())

    for i in 1:N_b
        b = b_grid[i] # Current state: cake stock
        
        # Define the action space (consumption c): [c_min_bound, b]
        N_c = 500 # High number for a precise consumption grid
        
        # Consumption can go from the lower bound up to the current stock b
        c_min_action = min(b, c_min_bound)
        c_grid = LinRange(c_min_action, b, N_c)
        
        # Handle tiny b values
        if b < c_min_bound
            c_grid = [b] 
        end
        
        max_val = -Inf
        best_c = 0.0

        # Maximization step: Loop over consumption choices
        for c in c_grid
            b_next = b - c # Next period's cake stock
            
            # Evaluate V(b') using the interpolation of V_old
            V_next_period = V_old_interp(b_next)
            
            # Bellman equation: u(c) + beta * V(b')
            current_val = u(c, c_min_bound) + beta * V_next_period

            # Update maximum value and optimal consumption
            if current_val > max_val
                max_val = current_val
                best_c = c
            end
        end

        # Store results
        V_new[i] = max_val
        policy[i] = best_c
    end

    # Check convergence
    max_diff = maximum(abs.(V_new - V_old))
    V_old .= V_new # Update V_old
    iter += 1
end

println("VFI converged in $iter iterations. Max difference: $max_diff")

# --- 4. Plot the Policy Function c(b) ---

# Get the analytical policy for comparison (c(b) = 0.1 * b)
analytical_policy = 0.1 * b_grid

policy_plot = plot(b_grid, policy,
                   label = "Computational Policy (VFI)",
                   xlabel = "Cake Stock (b)",
                   ylabel = "Consumption (c)",
                   title = "Optimal Consumption Policy (VFI vs. Analytical)",
                   legend = :topleft)

plot!(policy_plot, b_grid, analytical_policy, 
      label = "Analytical Policy (0.1b)", 
      linestyle = :dash, 
      linecolor = :red)

