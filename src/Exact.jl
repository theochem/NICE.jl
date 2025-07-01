using LinearAlgebra
using SciMLBase
using ForwardDiff
using NonlinearSolve

"""
    objective_function(f, rxn_system, K_eqs, xi)

Computes the objective function f_r(ξ) for each reaction r using the equation:
f_r = K^(r) * ∏(a_i^(-ν_i^(r))) - ∏(a_i^(ν_i^(r))) = 0
where the first product is over reactants (ν_i^(r) < 0) and the second over products (ν_i^(r) > 0).
Returns a vector of length R, where f_r = 0 at equilibrium.

Arguments:
- f::Vector{Float64}: The output vector to store the computed objective function values.
- rxn_system::ReactionSystem: The reaction system with stoichiometry, initial concentrations, and equilibrium constants.
- K_eqs::AbstractVector{Float64}: The equilibrium constants for each reaction, length R.
- xi::AbstractVector{T}: The reaction progress vector ξ, length R, where T is a subtype of Real (e.g., Float64 or Dual).
"""
function objective_function!(f::AbstractVector, rxn_system::ReactionSystem, K_eqs::AbstractVector{Float64}, xi::AbstractVector{T}) where {T<:Real}
    activities = rxn_system.concs_init .+ rxn_system.stoich * xi

    for r in 1:rxn_system.n_reaction
        term1 = K_eqs[r]
        term2 = one(T)

        for i in 1:rxn_system.n_species
            nu_ir = rxn_system.stoich[i, r]
            if nu_ir < 0
                term1 *= activities[i]^(-nu_ir)
            elseif nu_ir > 0
                term2 *= activities[i]^(nu_ir)
            end
        end
        f[r] = term1 - term2
    end
end
"""
    jacobian_function!(J, rxn_system, K_eqs, xi)

Computes the Jacobian matrix J of the nonlinear system in-place.

This implements the analytical Jacobian for the system of equations defined by
f_r = K^(r) * ∏(a_i^(-ν_i^(r))) - ∏(a_i^(ν_i^(r))) = 0
where the Jacobian element J_rs = ∂f_r/∂ξ_s is computed using:

J_rs = K^(r) * ∑_{i=1}^N [∏_{j=1,j≠i}^N (a_j^(0) + ∑_{r=1}^R ν_j^(r) ξ^(r))^(-ν_j^(r))] * [-ν_i^(r) * (a_i^(0) + ∑_{r=1}^R ν_i^(r) ξ^(r))^(-ν_i^(r)-1)] * (ν_i^(s))
    - ∑_{i=1}^N [∏_{j=1,j≠i}^N (a_j^(0) + ∑_{r=1}^R ν_j^(r) ξ^(r))^(ν_j^(r))] * [ν_i^(r) * (a_i^(0) + ∑_{r=1}^R ν_i^(r) ξ^(r))^(ν_i^(r)-1)] * (ν_i^(s))

Arguments:
- J::Matrix{Float64}: The output Jacobian matrix to be filled in-place.
- rxn_system::ReactionSystem: The reaction system containing stoichiometry and initial concentrations.
- K_eqs::AbstractVector{Float64}: The equilibrium constants for each reaction, length R.
- xi::Vector{Float64}: The reaction progress vector ξ, length R.


"""
function jacobian_function!(J::Matrix{Float64}, rxn_system::ReactionSystem, K_eqs::AbstractVector{Float64}, xi::Vector{Float64})
    R = rxn_system.n_reaction
    N = rxn_system.n_species

    activities = rxn_system.concs_init .+ rxn_system.stoich * xi

    if any(x -> x <= 0, activities)
        J .= Inf
        return nothing
    end

    P_reactants = ones(R)
    P_products = ones(R)
    for r in 1:R
        for i in 1:N
            nu_ir = rxn_system.stoich[i, r]
            if nu_ir < 0
                P_reactants[r] *= activities[i]^(-nu_ir)
            elseif nu_ir > 0
                P_products[r] *= activities[i]^(nu_ir)
            end
        end
    end

    for r in 1:R, s in 1:R
        sum_deriv1 = 0.0
        for i in 1:N
            if rxn_system.stoich[i, r] < 0
                sum_deriv1 += (-rxn_system.stoich[i, r] / activities[i]) * rxn_system.stoich[i, s]
            end
        end
        term1_deriv = K_eqs[r] * P_reactants[r] * sum_deriv1

        sum_deriv2 = 0.0
        for i in 1:N
            if rxn_system.stoich[i, r] > 0
                sum_deriv2 += (rxn_system.stoich[i, r] / activities[i]) * rxn_system.stoich[i, s]
            end
        end
        term2_deriv = P_products[r] * sum_deriv2

        J[r, s] = term1_deriv - term2_deriv
    end
end


"""
    solve(rxn_system::ReactionSystem; maxiters, abstol, reltol)

Solves for the equilibrium state of the given `ReactionSystem`.

This function sets up and solves the system of nonlinear equations for the
reaction extents (ξ) that bring all reactions to equilibrium. It uses the
`NonlinearSolve.jl` package with an analytically provided Jacobian for
efficiency and robustness.

Upon successful convergence, the `concs` field of the `rxn_system` object
is updated with the final equilibrium concentrations.

# Arguments
- `rxn_system::ReactionSystem`: The reaction system to solve. All necessary
  parameters (`stoich`, `concs_init`, `keq_vals`) are taken from this object.

# Keyword Arguments
- `maxiters::Integer=100`: Maximum number of iterations for the solver.
- `abstol::Real=1.0e-9`: Absolute tolerance for the solver.
- `reltol::Real=0.0`: Relative tolerance for the solver.

# Returns
- `solution`: The full solution object from `NonlinearSolve.solve`, which contains
  the final reaction extents (`solution.u`), the return code (`solution.retcode`),
  and other diagnostic information.
"""
function solve(
    rxn_system::ReactionSystem,
    K_eqs::AbstractVector{Float64};
    maxiters::Integer=1000,
    abstol::Real=1.0e-9,
    reltol::Real=0.0,
)
    # Define wrappers that unpack the parameter tuple p = (rxn_system, K_eqs)
    f_obj = (du, u, p) -> objective_function!(du, p[1], p[2], u)
    f_jac = (J, u, p) -> jacobian_function!(J, p[1], p[2], u)

    nls_func = NonlinearFunction(f_obj; jac=f_jac)
    u0 = zeros(Float64, rxn_system.n_reaction)
    params = (rxn_system, K_eqs)
    problem = NonlinearProblem(nls_func, u0, params)
    solution = NonlinearSolve.solve(problem, TrustRegion(); maxiters=maxiters, abstol=abstol, reltol=reltol)

    if SciMLBase.successful_retcode(solution)
        rxn_system.concs .= rxn_system.concs_init .+ rxn_system.stoich * solution.u
    end

    return solution
end

"""
    solve(rxn_system, K_eqs; max, φ=1.0, rate_consts=:forward, maxiters=1000, abstol=1.0e-9, reltol=0.0)

Solve the system deterministically.

This method solves a linear least squares problem to get an initial guess for
``\\mathbf{\\xi}``,
```math
{\\left(c_{mn}\\right)}^T \\left(\\xi_m^{\\left(0\\right)}\\right) =
    \\left(a_n - a_n^{\\left(0\\right)}\\right)
```
and then optimizes the equation for the simultaneous equilibria,
```math
K_\\text{eq}^{\\left(m\\right)} =
    \\prod_{n=1}^N {\\left( a_n^{\\left(0\\right)} +
    \\sum_{m=1}^M c_{mn} \\xi_m \\right)}^{c_{mn}}
```
using Newton's method + Trust Region. Instead of specifying the ``K_\\text{eq}``
values directly, they are approximated here using the forward or reverse rate
constants and a parameter ``\\phi``. If `:forward` is chosen, equilibrium constants
are derived as ``K_\\text{eq} = (k_f - 2\\phi)/(k_f - \\phi)``; if `:reverse`, as
``K_\\text{eq} = (k_r + \\phi)/k_r``. The function updates `rxn_system.concs`
upon successful convergence and returns a solution object.
"""
function solve(
    rxn_system::ReactionSystem;
    φ::Real=1.0,
    rate_consts::Symbol=:forward,
    maxiters::Integer=1000,
    abstol::Real=1.0e-9,
    reltol::Real=0.0,
)
    if rate_consts === :forward
        # Evaluate K_eqs from forward rate constants and φ_fwd
        K_eqs = (
            (rxn_system.fwd_rate_consts .- 2 * φ) ./
            (rxn_system.fwd_rate_consts - φ)
        )
    elseif rate_consts === :reverse
        # Evaluate K_eqs from reverse rate constants and φ_rev
        K_eqs = (rxn_system.rev_rate_consts .+ φ) ./ rxn_system.rev_rate_consts
    else
        throw(ArgumentError("invalid rate_consts value $rate_consts \
                             (must be :forward or :reverse)"))
    end
    solve(rxn_system, K_eqs; maxiters=maxiters, abstol=abstol, reltol=reltol)
end

"""
    hybrid_solve(rxn_system, K_eqs; n_iter=Int(1e+8), chunk_iter=Int(1e+4), ε=1.0e-4, ε_scale=1.0, ε_concs=0.0, tol_ε=0.0, maxiters=1000, abstol=1.0e-9, reltol=0.0)

Combines stochastic simulation and deterministic solving to find equilibrium
concentrations. First, it runs a Net-Event Kinetic Monte Carlo simulation via
`simulate` to approach equilibrium, then refines the result using the
deterministic `solve` method with provided ``K_\\text{eq}`` values. The function
updates `rxn_system.concs` and returns the final nonlinear solution object.
"""
function hybrid_solve(
    rxn_system::ReactionSystem,
    K_eqs::AbstractVector{Float64};
    n_iter::Integer=Int(1e+8),
    chunk_iter::Integer=Int(1e+4),
    ε::Real=1.0e-4,
    ε_scale::Real=1.0,
    ε_concs::Real=0.0,
    tol_ε::Real=0.0,
    maxiters::Integer=1000,
    abstol::Real=1.0e-9,
    reltol::Real=0.0,
)
    simulate(rxn_system; n_iter=n_iter, chunk_iter=chunk_iter, ε=ε, ε_scale=ε_scale, ε_concs=ε_concs, tol_ε=tol_ε)
    solve(rxn_system, K_eqs; maxiters=maxiters, abstol=abstol, reltol=reltol)
end