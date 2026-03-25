# Part 1, Exercise 2: Parameters

# ex2_parameters.jl
using Parameters

Params = @with_kw (
    nb = 50,
    b_min = -0.5,
    b_max = 5.0,
    yL = 1.0,
    yH = 2.0,
    β = 0.85,
    γ = 1.5,
    r = 0.10,
    tol_vfi = 1e-6,
    tol_kfe = 1e-6,
    max_iter = 500
)

trans_matrix(ρ) = [ρ (1-ρ); (1-ρ) ρ]