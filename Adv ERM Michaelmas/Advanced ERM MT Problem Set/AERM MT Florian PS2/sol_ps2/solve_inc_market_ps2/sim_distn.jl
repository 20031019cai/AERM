
# build a function that solves for the stationary distribution of the model

function solve_g(;sol_v=sol_v,x_int=x_int,params=params)
    @unpack tol_KFE,max_iter_KFE = params
    
    # unpack solution value function 
    c_sol        = sol_v.c_sol
    v_vec        = sol_v.v_vec
    B_next       = sol_v.B_next
    B_num_all = size(B_next)[1] # note B_num_all includes different y states
    B_grid_large = sol_v.B_grid_large
    y_grid_large = sol_v.y_grid_large
    pi_mat       = sol_v.pi_mat
    y_num = size(pi_mat)[1]
    
    # set up a large transition matrix 

    # the probability of going from B to B' is zero except for the case where B' = B_next
    trans_mat = zeros(B_num_all,B_num_all)
    for i = 1:B_num_all
        for j = 1:B_num_all
            if B_next[i] == B_grid_large[j]
                trans_mat[i,j]=1
            end
        end
    end

    # now add the stochastic income process
    B_num_small = Int(B_num_all/2)
    y_mat = kron(pi_mat,ones(B_num_small,B_num_small))
    trans_mat = trans_mat.*y_mat

    # sanity checks are always a good idea 
    findmax(sum(trans_mat, dims = 2)) # should be one
    findmin(sum(trans_mat, dims = 2)) # should be one

    # guess an initial distribution
    if x_int == 0 
         x_int = ones(size(sol_v.B_grid_large)[1])./size(sol_v.B_grid_large)[1] # uniform cause why not?
    else 
        x_int = x_int
    end
 
    x = copy(x_int)
    for k = 1:max_iter_KFE
        x_old = copy(x)
        x = trans_mat' * x
        if findmax(abs.(x-x_old))[1] < tol_KFE
            println("CONV stationary distn. after $k iterations.")    
            break
        end
    end

    
    return (trans_mat=trans_mat, x_ss=x)
end

