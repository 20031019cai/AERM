using(LinearAlgebra) # almost always relevant (but actually no needed here)
using(Parameters)    # for @with_kw, nice to keep track of an change parameters
using(Plots)         # for plotting
using(TickTock)      # for timing the code

# load a predifined function 
include("solve_v.jl") # solves the value function
include("def_params_baseline.jl") # solves the value function
include("sim_distn.jl") # solves the value function
include("comp_K_S_D.jl") # solves demand and supply for capital

tick() 
sol_v = solve_v(params=params_bl(),v_int=0) # comes from the code in "solve_v.jl"
tock()
tick()
sol_v = solve_v(params=params_bl(),v_int=sol_v.v_vec)
tock()

# plot the consumption function 
c_bl = plot(sol_v.c_sol, label = "c")

savefig("figures/c_a_bl.png")   # Export as PNG

# version with more grid points and larger upper bound
params_largegrid = params_bl(B_num=500)
tick()
sol_v_lg = solve_v(params=params_largegrid,v_int=0) # comes from the code in "solve_v.jl"
tock()
c_lg = plot(sol_v_lg.c_sol, label = "c")
savefig("figures/c_a_lg.png")   # Export as PNG

# version with more grid points and larger upper bound

params_largegrid_Bup = params_bl(B_num=1000,B_upper=10)
tick()
sol_v_lg_bup = solve_v(params=params_largegrid_Bup,v_int=0) # comes from the code in "solve_v.jl"
tock()
c_lg = plot(sol_v_lg_bup.c_sol, label = "c")
savefig("figures/c_a_lg_bup.png")   # Export as PNG
plot(c_lg)
# solve for the stationary distribution 
(mat,x) = solve_g(;sol_v=sol_v_lg_bup,x_int=0,params=params_largegrid_Bup) 



#sanity checks are always a good idea 
findmax(sum(mat, dims = 2)) # should be one
findmin(sum(mat, dims = 2)) # should be one


# compute aggregate demand for assets 
sol_K = comp_K_S_D(;dist_ss=x,params = params_largegrid_Bup) # the semicolon is important here, it makes sure that the function can only be called with named arguments

# doesn't match at all -- increase r 
r_low = params_largegrid_Bup.r



# doesn't match at all -- decrease r
r_ss=[];K_ss=[]
r_up = .10 
r_low = .08
B_num_final_sim = 2000
B_upper_final_sim = 16
params_next = params_bl(B_num=B_num_final_sim,B_upper=B_upper_final_sim,r=(r_up+r_low)/2)
sol_v = solve_v(params=params_next,v_int=0) # comes from the code in "solve_v.jl"
(mat,x) = solve_g(;sol_v=sol_v,x_int=0,params=params_next) 

# comment: here using solutions from the previous iteration as initial guess for the next iteration 
#          is key to make this routine faster. 

tick()
for k = 1:100
    sol_v = solve_v(params=params_next,v_int=sol_v.v_vec) # comes from the code in "solve_v.jl"
    (mat,x) = solve_g(;sol_v=sol_v,x_int=x,params=params_next) 
    sol_K = comp_K_S_D(;dist_ss=x,params = params_next)
    
    if sol_K.K_D > sol_K.K_S
        r_low = (r_up+r_low)/2    
        r_up = r_up
        params_next = params_bl(B_num=B_num_final_sim,B_upper=B_upper_final_sim,r=(r_low+r_up)/2)
    else
        r_low = r_low    
        r_up = (r_up+r_low)/2
        params_next = params_bl(B_num=B_num_final_sim,B_upper=B_upper_final_sim,r=(r_low+r_up)/2)
    end
    println("r_low = ",r_low,", r_up = ",r_up)
    println("sol_K.K_D-sol_K.K_S = ",sol_K.K_D-sol_K.K_S)
    if abs.(sol_K.K_D-sol_K.K_S)<params_next.tol_r
        r_ss = (r_low+r_up)/2
        K_ss = (sol_K.K_D+sol_K.K_S)/2
        println("Conv K,r.")
        break 
    end
end
tock()

# plot stationary distribution, equilibrium interest rate, and aggregate capital
c_y_b_final = plot(sol_v.B_grid_large[1:sol_v.B_num],sol_v.c_sol[1:sol_v.B_num],label="y_low",xlabel="Assets", ylabel="Consumption")
              plot!(sol_v.B_grid_large[sol_v.B_num+1:end],sol_v.c_sol[sol_v.B_num+1:end],label="y_high")

savefig("figures/c_y_b_final.png")   # Export as PNG
plot(c_y_b_final)


b_dist_ss = plot(sol_v.B_grid_large[1:sol_v.B_num],x[1:sol_v.B_num],label="y_low",xlabel="Assets", ylabel="pmf")
            plot!(sol_v.B_grid_large[1+sol_v.B_num:end],x[1+sol_v.B_num:end],label="y_high")

savefig("figures/b_dist_ss_final.png")   # Export as PNG
plot(b_dist_ss)

## big speed gains if using previous solution to v as initial guess!
println("*** end of code ***")






