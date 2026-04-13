# Part 1, Exercise 2: Value Function Iteration

# ex2_vfi.jl

function solve_household(p, ρ::Float64)
    
    # Unpack

        nb = p.nb
        b_grid = collect(LinRange(p.b_min, p.b_max, nb))
        y_grid = [p.yL, p.yH]
        ny = 2
        Π = trans_matrix(ρ)

    # CRRA utility

        function u(c)
            c > 0 ? (c^(1 - p.γ)) / (1 - p.γ) : -1e10
        end

    # Pre-compute consumption for all (b, y, b') triples
    # c(b, y, b') = (1+r)*b + y - b'
    # Store as matrix: rows = (b,y) states, cols = b' choices

        n_states = nb * ny
        c_mat = zeros(n_states, nb)
        for iy in 1:ny
            for ib in 1:nb
                idx = (iy - 1) * nb + ib
                for ib_next in 1:nb
                    c_mat[idx, ib_next] = (1 + p.r) * b_grid[ib] + y_grid[iy] - b_grid[ib_next]
                end
            end
        end

    # Utility matrix

        u_mat = u.(c_mat)

    # Initialise value function

        V = -ones(n_states)
        V_new = similar(V)

        # stores index of optimal b'
        policy_idx = zeros(Int, n_states) 

    # VFI loop

    converged = false
    for iter in 1:p.max_iter
        # Compute expected continuation value: EV(b') = Σ_y' Π(y,y') V(b',y')
        # V is ordered as [V(b1,y1), V(b2,y1),...,V(bN,y1), V(b1,y2),...,V(bN,y2)]
        
        # nb × ny
            V_mat = reshape(V, nb, ny)          
        
        # nb × ny (expected V for each b' and current y)
            EV_mat = V_mat * Π'                 
        
        # back to vector, same ordering as V
            EV = vec(EV_mat)                    

        # Maximise over b' for each state (b, y)
        
        for iy in 1:ny

            # EV(b', y) for current y
                EV_slice = EV[(iy-1)*nb+1 : iy*nb]  
                for ib in 1:nb
                    idx = (iy - 1) * nb + ib
                    best_val = -Inf
                    best_ib = 1
                    for ib_next in 1:nb
                        val = u_mat[idx, ib_next] + p.β * EV_slice[ib_next]
                        if val > best_val
                            best_val = val
                            best_ib = ib_next
                        end
                    end
                    V_new[idx] = best_val
                    policy_idx[idx] = best_ib
                end
        end

    # Check convergence
        
        if maximum(abs.(V_new - V)) < p.tol_vfi
            V .= V_new
            converged = true
            println("VFI converged after $iter iterations (ρ=$ρ)")
            break
        end
        V .= V_new
    end

    if !converged
        println("VFI did NOT converge for ρ=$ρ")
    end

    # Extract policy functions
        b_next = [b_grid[policy_idx[i]] for i in 1:n_states]
        c_policy = [(1 + p.r) * b_grid[mod1(i, nb)] + y_grid[div(i-1, nb)+1] - b_next[i] for i in 1:n_states]

    return (V=V, policy_idx=policy_idx, b_next=b_next, c_policy=c_policy,
            b_grid=b_grid, y_grid=y_grid, nb=nb, ny=ny, converged=converged)
end