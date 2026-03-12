# run phelps model of entrepreneurial risk ß

using(LinearAlgebra) # almost always relevant (but actually no needed here)
using(Parameters)    # for @with_kw, nice to keep track of an change parameters
using(Plots)         # for plotting
using(TickTock)      # for timing the code

# load a predefined function 
include("make_discrete_asset_returns.jl") # solves the value function
include("solve_v.jl") # solves the value function
include("def_params_baseline.jl") # solves the value function



params_bl = parameters_bl()
(R_vec,trans_mat) = gen_R_trans_mat(;params = params_bl)

sol_v_bl = solve_v(;params = params_bl,v_int=0,trans_mat=trans_mat,R_vec=R_vec)
# plot consumption policy function for high and low interest rate, relative to normalized wealth
c_x = plot(sol_v_bl.c_sol[1:params_bl.x_num],label="R low")
plot!(sol_v_bl.c_sol[1+params_bl.x_num:end],label="R high")
savefig("figures/c_x.png")   # Export as PNG
# plot policy function for high and low interest rate relative to wealth 
c_sol = reshape(sol_v_bl.c_sol,params_bl.x_num,2)
w_grid = reshape(sol_v_bl.w,params_bl.x_num,2)
v_vec = reshape(sol_v_bl.v_vec,params_bl.x_num,2)

c_w = plot(w_grid[:,1],c_sol[:,1],lw=3, linecolor=:red , linestyle= :dash,label="R low")
plot!(w_grid[:,2],c_sol[:,2],lw=1, linecolor=:blue, label = "R high")
plot!(w_grid[:,2],w_grid[:,2]*(1-params.β),lc=:black,ls=:dot,label="w(1-β)")
savefig("figures/c_w.png")


###########################
# explore with iid returns
###########################

params_iid = parameters_bl(ρ_L = .5,  ρ_H = .5)
(R_vec,trans_mat) = gen_R_trans_mat(;params = params_iid)
sol_v_iid = solve_v(;params = params_iid,v_int=0,trans_mat=trans_mat,R_vec=R_vec)

# plot consumption policy function for high and low interest rate, relative to normalized wealth
c_iid = plot(sol_v_iid.c_sol[1:params_iid.x_num],label="R low")
plot!(sol_v_iid.c_sol[1+params_iid.x_num:end],label="R high")
savefig("figures/c_x_iid.png")   # Export as PNG

# plot policy function for high and low interest rate relative to wealth 
c_sol = reshape(sol_v_iid.c_sol,params_iid.x_num,2)
w_grid = reshape(sol_v_iid.w,params_iid.x_num,2)
v_vec = reshape(sol_v_iid.v_vec,params_iid.x_num,2)

c_w_iid = plot(w_grid[:,1],c_sol[:,1],lw=3, linecolor=:red , linestyle= :dash,label="R low")
plot!(w_grid[:,2],c_sol[:,2],lw=1, linecolor=:blue, label = "R high")
plot!(w_grid[:,2],w_grid[:,2]*(1-params.β),lc=:black,ls=:dot,label="w(1-β)")
savefig("figures/c_w_iid.png")


## compute theoretical value function
function comp_v_theory(w,β)
    v = log.(w)./(1 - β) .+ log(1-β)/(1-β) .+ (β/(1-β)^2)*(log(β)) 
    # this is true because expected log return zero
    return v 
end 

vv = comp_v_theory(w_grid[:,2],params_iid.β)
v_iid = plot(vv,v_vec[:,2],label = "",lw=2,lc=:blue, ylabel="v-iteration", xlabel="v-exact")
plot!(vv,vv,lc = :red, ls =:dash, label="45 degree line") # mimicks 45 degree line
savefig("figures/v_iid.png")


