using(LinearAlgebra) # almost always relevant (but actually no needed here)
using(Parameters)    # for @with_kw, nice to keep track of an change parameters
using(Plots)         # for plotting
using(TickTock)      # for timing the code

# load a predifined function 
include("solve_v.jl") # solves the value function
include("def_params_baseline.jl") # solves the value function
include("sim_distn.jl") # solves the value function
include("comp_K_S_D.jl") # solves demand and supply for capital


################################
# 1. Stationary income process # 
################################

params_ex1 = params_bl()
params_ex1_highp = params_bl(p_y = 0.99)  
pi_mat = [params_ex1.p_y (1-params_ex1.p_y); 1-params_ex1.p_y params_ex1.p_y]
pi_mat_hp = [params_ex1_highp.p_y (1-params_ex1_highp.p_y); (1-params_ex1_highp.p_y) params_ex1_highp.p_y]
 
x_int = [1;0]
year =1
x_mat = x_int'
while year <= 50
    x_mat = [x_mat; (pi_mat'*x_mat[year,:])']
    year +=1
end
x_int = [1;0]
year =1
x_mat_hp = x_int'
while year <= 50
    x_mat_hp = [x_mat_hp; (pi_mat_hp'*x_mat_hp[year,:])']
    year +=1
end
plot_share_low  = plot([0:1:50],x_mat[:,1],label="share low income",xlabel="year")
plot!([0:1:50],x_mat_hp[:,1],label="share low income with high persistence",xlabel="year")
savefig("figures/plot_share_low.png")   # Export as PNG


##################################################################
# 2. solve household problem for different values of persistence #
##################################################################

# values over which I loop 
val_ρ = [.05,.25,.5,.6,.7,.8,.9,.95,.98,.99]
# solve value function
sol_v_vec = [solve_v(params = params_bl(p_y = ρ), v_int=0) for ρ in val_ρ]
# solve stationary distribution
solv_g_vec = [solve_g(;sol_v=sol_v,x_int=0,params=params_bl()) for sol_v in sol_v_vec]
B_agg_vec = zeros(size(val_ρ)[1])
for i =1:size(val_ρ)[1]
    B_agg_vec[i] = (solv_g_vec[i].x_ss)'*sol_v_vec[i].B_grid_large
    println((solv_g_vec[i].x_ss)'*sol_v_vec[i].y_grid_large)
end


# plot aggregate demand for assets
plot_asset_rho  = plot(val_ρ,B_agg_vec,label="asset demand",xlabel="rho")
savefig("figures/plot_asset_rho.png")   # Export as PNG


## changing borrowing constraint ##
val_B_min = collect(LinRange(-params_bl().y_low/params_bl().r,0,5))

# solve value function
sol_v_vec = [solve_v(params = params_bl(B_lower = Bmin), v_int=0) for Bmin in val_B_min]
# solve stationary distribution
solv_g_vec = [solve_g(;sol_v=sol_v,x_int=0,params=params_bl()) for sol_v in sol_v_vec]
B_agg_vec = zeros(size(val_B_min)[1])
for i =1:size(val_B_min)[1]
    B_agg_vec[i] = (solv_g_vec[i].x_ss)'*sol_v_vec[i].B_grid_large
    println((solv_g_vec[i].x_ss)'*sol_v_vec[i].y_grid_large)
end


# plot aggregate demand for assets
plot_asset_bconstraint  = plot(val_B_min,B_agg_vec,label="asset demand",xlabel="B_min")
savefig("figures/plot_asset_bmin.png")   # Export as PNG



############################################
# do the whole thing with endogenous wages #
############################################


# compute aggregate demand for assets 

# doesn't match at all -- increase r 

# doesn't match at all -- decrease r
r_ss=[];K_ss=[]
r_up = .095 
r_low = .091
B_num_final_sim = 1000
B_upper_final_sim = 13
params_next = params_bl(B_num=B_num_final_sim,B_upper=B_upper_final_sim,r=(r_up+r_low)/2)
sol_v = solve_v_wage(params=params_next,v_int=0) # comes from the code in "solve_v.jl"
(mat,x) = solve_g(;sol_v=sol_v,x_int=0,params=params_next) 

# comment: here using solutions from the previous iteration as initial guess for the next iteration 
#          is key to make this routine faster. 

tick()
for k = 1:100
    sol_v = solve_v_wage(params=params_next,v_int=sol_v.v_vec) # comes from the code in "solve_v.jl"
    (mat,x) = solve_g(;sol_v=sol_v,x_int=x,params=params_next) 
    sol_K = comp_K_S_D_wage(;dist_ss=x,params = params_next,w_t=sol_v.w_t)
    
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


println("*** end of code ***")






