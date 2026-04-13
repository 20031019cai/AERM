# # Usage example:
# npts = 10
# rho = 0.9
# sig = 0.01
# gridpts, transition = rouwenhorst(npts, rho, sig)



function gen_R_trans_mat(;params = params)
    @unpack R_high, R_low, ρ_L,ρ_H = params
    # define initial matrix
    mat = [ρ_L 1-ρ_L; (1-ρ_H) ρ_H]
  
    R_vec = [R_low;R_high]
    return (R_vec=R_vec,trans_mat=mat)

end