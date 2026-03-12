
# this file tests github autopilot 
function solve_v(;params = params,v_int=v_int) # the semicolon is important here, it makes sure that the function can only be called with named arguments
                                   # you can also do solve_9(50,0,4,.8,...) but you have to remember the order of the arguments, much can go wrong...         
    # unpack parameters
    @unpack B_num,B_lower,B_upper,y_low,y_high,p_y,β,r,γ,max_iter_HJB,tol_HJB = params

    # set up grid space
    B_grid = collect(LinRange(B_lower,B_upper,B_num))
    y_grid = [y_low;y_high]
    pi_mat = [p_y 1-p_y;1-p_y p_y]
    y_num = size(y_grid)[1]

    convergence_v = false # boolean to check if converged
    
    # initialize value function -- could be any finite number, and doesn't have to be the same
    # if you want to converge quick, a good guess helps
    if v_int == 0
        v_0 = -ones(y_num*B_num) # states_y X states_B
    else 
        v_0 = v_int
    end
    v_vec = copy(v_0)
    E_v = reshape(reshape(v_vec,B_num,y_num)*pi_mat',B_num*y_num)
    policy_function = zeros(B_num*y_num) # this one will be an index, where the index refers to the optimal choice of B_next
    
    # make payoff vector -- will be of size B_num*B_num*y_num (note ``curse of dimensionality'') 
    # quick function to deal with negative consumption, important note: γ is defined above, and then used inside this function
    function u(c)
        if c>0
            return (c^(1-γ))/(1-γ)
        else
            return -1000000  # we just want something very negative here, so that in the end no one is ever gonna pick a negative consumption level
        end
    end

    y_grid_large = kron(y_grid,ones(B_num*B_num))
    c_vec = kron(ones(y_num),(1+r)*kron(B_grid,ones(B_num))-kron(ones(B_num),B_grid)) + y_grid_large # when prepping, I made a mistake here and it took me two fours to find the bug because I didn't thought it possible to make a simple mistake so early on
    u_c_vec = u.(c_vec)
    
    ## start the iterative procedure -- maximum iter is max_iter_HJB
    for l = 1:max_iter_HJB
        # blow up the value function from a B_num*y_num vector to a B_num*B_num*y_num vector, so you can take the max over all choices of B_next within specific B_current
        v_continue=[]
        for i=1:y_num
            v_continue = [v_continue; kron(ones(B_num),E_v[B_num*(i-1)+1:B_num*(i)])]
        end
        
        # combine the payoff vector with the continuation value (payoff always the same but v is being updated)
        v_vec_big = u_c_vec + β*v_continue
        v_vec_old = copy(v_vec) # copy matters here instead of simply setting v_vec_old = (v_vec), figure out why
 
        for i=1:B_num*y_num # we have to iterate over all states (B*y)
            v_vec[i] = findmax(v_vec_big[1+B_num*(i-1):B_num*i])[1] # findmax returns a tuple (max_value, index), we want the first element of that tuple    
        end
        
        # get the expected value -- pi_mat is the stochastic matrix that matters now
        E_v = reshape(reshape(v_vec,B_num,y_num)*pi_mat',B_num*y_num) # you might be able to write this matrix multiplication more elegantly, but make sure you multiply v the right way 
        if findmax(abs.(v_vec-v_vec_old))[1]< tol_HJB # check for convergence
            # if converged, get the policy function and break the loop
            for i=1:B_num*y_num 
                policy_function[i] = Int(findmax(v_vec_big[1+B_num*(i-1):B_num*i])[2])
            end 
            convergence_v = true # boolean to check if converged

            break
        end
    end 
    # turn policy function from index into actual consumption
#    policy_function=Int.(policy_function)
#    B_next = []
#    for i=1:B_num*y_num
#        B_next = [B_next; B_grid[policy_function[i]]]
#    end
println("convergence_v: ",convergence_v)
    B_next = [B_grid[Int(policy_function[i])] for i in 1:(B_num * y_num)]
    c_sol = (1+r)*kron(ones(y_num),B_grid)-B_next + kron(y_grid,ones(B_num))
    B_grid_large = kron(ones(y_num),B_grid) # this is just for plotting purposes, and to remind us how we orderd the states encoded in the vector
    return (c_sol=c_sol,convergence_v=convergence_v,v_vec=v_vec,B_grid_large = B_grid_large,B_num=B_num,y_grid_large = kron(y_grid,ones(B_num)),pi_mat=pi_mat,B_next=B_next)
end