using LinearAlgebra

using Reexport

using SciMLBase
using ForwardDiff
using NonlinearSolve

using Evolutionary
@reexport using Evolutionary: ES, CMAES, GA, DE, NSGA2, TreeGP

export solve!
export evolutionary_solve!

struct ConvergenceError <: Exception
    msg::String
end

Base.showerror(io::IO, e::ConvergenceError) = print(io, e.msg)

"""
    solve!(rxn_system, K_eqs; maxiters=1000, abstol=1.0e-9, reltol=0.0)

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
using user-provided ``K_\\text{eq}`` values and the Newton's method + Trust
Region method.
"""
function solve!(
    rxn_system::ReactionSystem,
    K_eqs::Vector{Float64};
    maxiter::Int=1000,
    abstol::Float64=1.0e-9,
    reltol::Float64=0.0,
)
    # Define nonlinear problem (derives Jacobian via autodiff)
    problem = NonlinearSolve.NonlinearProblem(
        # Objective function
        (f, ξ, _) -> f_K_eqs(f, ξ, rxn_system, K_eqs),
        # Initial guess for reaction extents ξ
        rxn_system.stoich \ (rxn_system.concs - rxn_system.concs_init);
        # Solver options
        maxiters=maxiter,
        abstol=abstol,
        reltol=reltol,
    )
    # Run the nonlinear optimization
    solution = NonlinearSolve.solve(problem, NonlinearSolve.TrustRegion())
    # If nonlinear optimization was successful, update `rxn_system.concs`;
    # otherwise throw ConvergenceError exception
    if !SciMLBase.successful_retcode(solution)
        throw(ConvergenceError("Did not converge"))
    end
    rxn_system.concs .= rxn_system.concs_init
    rxn_system.concs .+= rxn_system.stoich * solution.u
    # Return the nonlinear solution
    solution
end

"""
    solve!(rxn_system, K_eqs; maxiters=1000, abstol=1.0e-9, reltol=0.0)

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
using Newton's method + Trust Region.

Instead of specifying the ``K_\\text{eq}`` values directly, they are
approximated here using the forward or reverse rate constants and a parameter
``\\phi``.
"""
function solve!(
    rxn_system::ReactionSystem;
    φ::Float64=1.0,
    rate_consts::Symbol=:forward,
    maxiter::Int=1000,
    abstol::Float64=1.0e-9,
    reltol::Float64=0.0,
)
    solve!(
        rxn_system,
        approximate_K_eqs(rxn_system, φ, rate_consts);
        maxiter=maxiter,
        abstol=abstol,
        reltol=reltol,
    )
end

"""
    evolutionary_solve!(rxn_system, K_eqs, algorithm; abstol=Inf, reltol=Inf, maxiter=1000)

Solve the system using an evolutionary algorithm.

Solves the nonlinear problem described by [`solve!`](@ref) by using one of the
evolutionary algorithms from
[Evolutionary.jl](https://wildart.github.io/Evolutionary.jl/stable/).
These algorithm types are reexported by NICE.jl.
"""
function evolutionary_solve!(
    rxn_system::ReactionSystem,
    K_eqs::Vector{Float64},
    algorithm::TOpt;
    abstol::Float64=Inf,
    reltol::Float64=Inf,
    maxiter::Int=1000,
) where {TOpt <: Evolutionary.AbstractOptimizer}
    f = similar(K_eqs)
    solution = Evolutionary.optimize(
        # Objective function
        ξ -> f_K_eqs_cma(f, ξ, rxn_system, K_eqs),
        # Initial guess for reaction extents ξ
        rxn_system.stoich \ (rxn_system.concs - rxn_system.concs_init),
        # Evolutionary algorithm
        algorithm,
        # Solver options
        Evolutionary.Options(iterations=maxiter, abstol=abstol, reltol=reltol),
    )
    # If nonlinear optimization was successful, update `rxn_system.concs`
    # otherwise throw ConvergenceError exception
    if !solution.converged
        throw(ConvergenceError("Did not converge"))
    end
    rxn_system.concs .= rxn_system.concs_init
    rxn_system.concs .+= rxn_system.stoich * solution.minimizer
    # Return the nonlinear solution
    solution
end

"""
    evolutionary_solve!(rxn_system, K_eqs, algorithm; abstol=Inf, reltol=Inf, maxiter=1000)

Solve the system using an evolutionary algorithm.

Solves the nonlinear problem described by [`solve!`](@ref) by using one of the
evolutionary algorithms from
[Evolutionary.jl](https://wildart.github.io/Evolutionary.jl/stable/).
These algorithm types are reexported by NICE.jl.

Instead of specifying the ``K_\\text{eq}`` values directly, they are
approximated here using the forward or reverse rate constants and a parameter
``\\phi``.
"""
function evolutionary_solve!(
    rxn_system::ReactionSystem,
    algorithm::TOpt;
    φ::Float64=1.0,
    rate_consts::Symbol=:forward,
    abstol::Float64=Inf,
    reltol::Float64=Inf,
    maxiter::Int=1000,
) where {TOpt <: Evolutionary.AbstractOptimizer}
    evolutionary_solve!(
        rxn_system,
        approximate_K_eqs(rxn_system, φ, rate_consts),
        algorithm;
        abstol=abstol,
        reltol=reltol,
        maxiter=maxiter,
    )
end

function approximate_K_eqs(
    rxn_system::ReactionSystem,
    φ::Float64,
    rate_consts::Symbol,
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
    K_eqs
end

function f_K_eqs(f, ξ, rxn_system, K_eqs)::Nothing
    for i in 1 : rxn_system.n_reaction
        f_i1 = K_eqs[i]
        f_i2 = 1.0
        for j in 1 : rxn_system.n_species
            # Reactants
            if rxn_system.stoich[j, i] < 0
                f_i1 *= (
                    rxn_system.concs_init[j] + rxn_system.stoich[j, :] · ξ
                ) ^ -rxn_system.stoich[j, i]
            # Products
            elseif rxn_system.stoich[j, i] > 0
                f_i2 *= (
                    rxn_system.concs_init[j] + rxn_system.stoich[j, :] · ξ
                ) ^ rxn_system.stoich[j, i]
            end
        end
        f[i] = f_i1 - f_i2
    end
    nothing
end

function f_K_eqs_cma(f::Vector{Float64}, ξ::Vector{Float64}, rxn_system::ReactionSystem, K_eqs::Vector{Float64})::Float64
    f_K_eqs(f, ξ, rxn_system, K_eqs)
    √(f · f)
end
