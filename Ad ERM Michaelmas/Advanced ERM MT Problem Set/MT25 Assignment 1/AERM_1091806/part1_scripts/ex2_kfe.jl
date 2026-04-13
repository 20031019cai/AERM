# Part 1, Exercise 2: Kolmogorov-Forward Equation
# kfe.jl — Compute stationary distribution over (b, y) states

function stationary_dist(sol, ρ::Float64, p)
    nb = sol.nb
    ny = sol.ny
    n_states = nb * ny
    Π = trans_matrix(ρ)

    # Build full transition matrix over (b, y) → (b', y')
    # For state i = (b, y): policy says go to b' = b_grid[policy_idx[i]]
    # Then y transitions according to Π
    T = zeros(n_states, n_states)

    for iy in 1:ny
        for ib in 1:nb
            idx_from = (iy - 1) * nb + ib
            ib_next = sol.policy_idx[idx_from]
            for iy_next in 1:ny
                idx_to = (iy_next - 1) * nb + ib_next
                T[idx_from, idx_to] = Π[iy, iy_next]
            end
        end
    end

    # Iterate to find stationary distribution: μ' = T' * μ
    μ = ones(n_states) ./ n_states  # uniform initial guess

    for k in 1:p.max_iter
        μ_old = copy(μ)
        μ = T' * μ
        if maximum(abs.(μ - μ_old)) < p.tol_kfe
            println("KFE converged after $k iterations")
            break
        end
    end

    return μ
end

# Compute aggregate asset supply: B = Σ μ(b,y) * b
function agg_asset_supply(μ, sol)
    b_grid_large = vcat(sol.b_grid, sol.b_grid)  # repeat for each y
    return dot(μ, b_grid_large)
end