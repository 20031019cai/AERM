# define parameters for the model

params_bl = @with_kw ( 
    A_tfp = 1,
    B_num        = 1000,
    B_lower      = 0,
    B_upper      = 10,
    y_low        = .5,
    y_high       = 1.5,
    p_y          = .8 ,
    β            = .9,
    γ            =  2,
    r            = .0915,
    δ            = .05,
    α            = .4,
    max_iter_HJB = 400,
    max_iter_KFE = 400,
    tol_HJB      = 10^-7,
    tol_KFE      = 10^-8,
    tol_r        = 10^-4,)
