# define parameters for the model

parameters_bl = @with_kw ( 
    x_num        = 2000,
    x_lower      = 0.0,
    x_upper      = 15,
    c_gov        =0.0001,
    R_high       = 2,
    R_low        = .5,
    ρ_L            = .8,
    ρ_H            = .6,
    β            = .8,
    max_iter_HJB = 400,
    tol_HJB      = 10^-8,)
