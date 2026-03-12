
# compute demand and supply for capital
function comp_K_S_D(;dist_ss=dist_ss,params = params) # the semicolon is important here, it makes sure that the function can only be called with named arguments
                                   # you can also do solve_9(50,0,4,.8,...) but you have to remember the order of the arguments, much can go wrong...         
    # unpack parameters
    @unpack B_num,B_lower,B_upper,y_low,y_high,p_y,β,r,α,δ,γ,max_iter_HJB,tol_HJB = params

    # set up grid space
    B_grid = collect(LinRange(B_lower,B_upper,B_num))
    y_grid = [y_low;y_high]
    pi_mat = [p_y 1-p_y;1-p_y p_y]
    y_num = size(y_grid)[1]
    B_grid_large = kron(ones(y_num),B_grid) 
    y_grid_large = kron(y_grid,ones(B_num))
    # from household problem
    K_supply = dist_ss'* B_grid_large
    Lw = dist_ss'*y_grid_large # this one here has no economics in it and is just a function of the income process
    
    # from firm problem
    K_demand = (α/(1-α))*(Lw/(r+δ))
   
    return (K_S=K_supply,K_D=K_demand,Lw=Lw)
end