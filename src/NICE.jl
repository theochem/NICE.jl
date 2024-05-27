"""
NICE (N-species Ice Table) module -- provides simultaneous equilibria solvers:
- [`simulate`](@ref) Net-event kinetic Monte Carlo solver
- [`solve`](@ref) Nonlinear equations solver (Newton trust region method)

The [`ReactionSystem`](@ref) type is provided as an interface for these solvers.
"""
module NICE

using LinearAlgebra

using SciMLBase
using ForwardDiff
using NonlinearSolve

export ReactionSystem, simulate, solve

"""
Definition of the simultaneous equilibria reaction system. This type is the
input to the solver functions [`simulate`](@ref) and [`solve`](@ref), and is
updated in-place.

The reaction system consists of
- ``M`` Reactions involving ``N`` species;
- a matrix ``C = \\left(c_{mn}\\right)`` defining the stoichiometry of the
  system; rows correspond to species and columns correspond to reactions; each
  element is the stoichoimetric coefficient of the ``n``th species in the
  ``m``th reaction, and is positive if the species is a product, and negative
  if it is a reactant;
- the forward and reverse rate constants ``k_f`` and ``k_r`` for each reaction;
- initial concentrations or activities ``\\mathbf{a} = \\left(a_n\\right)``
  of each species;
and can be written (through some abuse of notation) as 
```math
\\begin{aligned}
\\sum_{\\left. n \\in \\left[1, N\\right] \\middle| c_{1n} < 0 \\right.} &\\left|c_{1n}\\right| I_n
&\\xrightleftharpoons[k_r^{\\left(1\\right)}]{k_f^{\\left(1\\right)}}
\\sum_{\\left. n \\in \\left[1, N\\right] \\middle| c_{1n} > 0 \\right.} &\\left|c_{1n}\\right| I_n \\\\
\\vdots & ~ & ~ & \\vdots \\\\
\\sum_{\\left. n \\in \\left[1, N\\right] \\middle| c_{mn} < 0 \\right.} &\\left|c_{mn}\\right| I_n
&\\xrightleftharpoons[k_r^{\\left(M\\right)}]{k_f^{\\left(M\\right)}}
\\sum_{\\left. n \\in \\left[1, N\\right] \\middle| c_{Mn} > 0 \\right.} &\\left|c_{Mn}\\right| I_n
\\end{aligned}
```
"""
mutable struct ReactionSystem
    n_species::Int
    n_reaction::Int
    stoich::Matrix{Float64}
    rev_rate_consts::Vector{Float64}
    fwd_rate_consts::Vector{Float64}
    concs_init::Vector{Float64}
    n_iter::Int
    time::Float64
    concs::Vector{Float64}
    rev_rates::Vector{Float64}
    fwd_rates::Vector{Float64}
    net_rates::Vector{Float64}
end

"""
    ReactionSystem(stoich, concs, rev_rate_consts, fwd_rate_consts)

Instantiate a [`ReactionSystem`](@ref) instance using the forward and reverse
rate constants for each reaction.
"""
function ReactionSystem(stoich, concs, rev_rate_consts, fwd_rate_consts)
    n_species = size(stoich, 1)
    n_reaction = size(stoich, 2)
    concs_init = copy(concs)
    n_iter = 0
    time = 0.
    concs = copy(concs)
    rev_rates = zeros(Float64, n_reaction)
    fwd_rates = zeros(Float64, n_reaction)
    net_rates = zeros(Float64, n_reaction)
    return ReactionSystem(
        n_species, n_reaction, stoich, rev_rate_consts, fwd_rate_consts, concs_init,
        n_iter, time, concs, rev_rates, fwd_rates, net_rates,
    )
end

"""
    ReactionSystem(stoich, concs, K_eqs; φ)

Instantiate a [`ReactionSystem`](@ref) instance using the equilibrium constants
``K_\\text{eq}^{\\left(m\\right)}`` to derive the forward and reverse rate
constants for each reaction:
```math
\\begin{aligned}
k_r &= \\frac{\\phi}{K_\\text{eq} - 1} \\\\
k_f &= \\phi - k_r
\\end{aligned}
```
"""
function ReactionSystem(stoich, concs, K_eqs; φ::Real=1.0)
    n_species = size(stoich, 1)
    n_reaction = size(stoich, 2)
    rev_rate_consts = φ ./ (K_eqs .+ 1)
    fwd_rate_consts = φ .- rev_rate_consts
    concs_init = copy(concs)
    n_iter = 0
    time = 0.
    concs = copy(concs)
    rev_rates = zeros(Float64, n_reaction)
    fwd_rates = zeros(Float64, n_reaction)
    net_rates = zeros(Float64, n_reaction)
    return ReactionSystem(
        n_species, n_reaction, stoich, rev_rate_consts, fwd_rate_consts, concs_init,
        n_iter, time, concs, rev_rates, fwd_rates, net_rates,
    )
end

"""
    simulate(rxn_system;
             n_iter=Int(1e+8), chunk_iter=Int(1e+4),
             ε=1.0e-4, ε_scale=1.0, ε_concs=0.0,
             tol_t=Inf, tol_ε=0.0, tol_concs=0.0)

Run a *N*et-*E*vent *K*inetic *M*onte *C*arlo (NEKMC) simulation to find the
equilibrium concentrations of the reaction system.
"""
function simulate(
    rxn_system;
    n_iter=Int(1e+8),
    chunk_iter=Int(1e+4),
    ε=1.0e-4,
    ε_scale=1.0,
    ε_concs=0.0,
    tol_t=Inf,
    tol_ε=0.0,
    tol_concs=0.0,
)
    # Allocate probability vector, Δconcs vector, Δt
    pvec = zeros(Float64, rxn_system.n_reaction)
    Δconcs = similar(rxn_system.concs)
    # Run simulation
    for _ in 1 : chunk_iter : n_iter
        Δtime = 0.0
        for _ in 1 : chunk_iter
            # Set Δconcs to concentration at step (n_iter - 1)
            Δconcs .= rxn_system.concs
            # Update rates of reaction
            update_rates(rxn_system)
            # Select and carry out random reaction
            i_rxn = select_reaction(rxn_system, pvec)
            Δtime += do_reaction(rxn_system, i_rxn, ε)
        end
        rxn_system.time += Δtime
        # Check for time convergence
        if Δtime > tol_t
            return :TimeConvergence
        end
        # Compute actual Δconcs at step (n_iter)
        Δconcs .-= rxn_system.concs
        norm_Δconcs = √(Δconcs · Δconcs)
        # Check for concentration convergence
        if norm_Δconcs < tol_concs
            return :ConcentrationConvergence
        # Check for decrease in concentration step size
        elseif norm_Δconcs < ε * ε_concs
            ε *= ε_scale
        end
        # Check for concentration step size convergence
        if ε < tol_ε
            return :StepSizeConvergence
        end
        rxn_system.n_iter += 1
    end
    :IterationLimit
end

function update_rates(rxn_system)
    for i in 1 : rxn_system.n_reaction
        rev_rate = rxn_system.rev_rate_consts[i]
        fwd_rate = rxn_system.fwd_rate_consts[i]
        for j in 1 : rxn_system.n_species
            s = rxn_system.stoich[j, i]
            if s >= 0.
                rev_rate *= rxn_system.concs[j] ^ s
            else
                # fwd_rate *= rxn_system.concs[j] ^ abs(s)
                fwd_rate *= rxn_system.concs[j] ^ (-s)
            end
        end
        rxn_system.rev_rates[i] = rev_rate
        rxn_system.fwd_rates[i] = fwd_rate
        rxn_system.net_rates[i] = fwd_rate - rev_rate
    end
end

function select_reaction(rxn_system, pvec)
    # Update probability vector
    p = 0.
    for i in 1 : rxn_system.n_reaction
        p += abs(rxn_system.net_rates[i])
        pvec[i] = p
    end
    # Select random reaction
    p *= rand(Float64)
    for i in 1 : rxn_system.n_reaction
        if pvec[i] > p
            return i
        end
    end
    # sentinel return value
    rxn_system.n_reaction
end

function do_reaction(rxn_system, i_rxn, ε)
    rate = rxn_system.net_rates[i_rxn]
    if rate >= 0
        for j = 1 : rxn_system.n_species
            rxn_system.concs[j] += rxn_system.stoich[j, i_rxn] * ε
        end
        return ε / rate
    else
        for j = 1 : rxn_system.n_species
            rxn_system.concs[j] -= rxn_system.stoich[j, i_rxn] * ε
        end
        return -ε / rate
    end
end

"""
    solve(rxn_system, K_eqs; maxiters=1000, abstol=1.0e-9, reltol=0.0)

Solve the system deterministically.
"""
function solve(rxn_system, K_eqs; maxiters=1000, abstol=1.0e-9, reltol=0.0)
    # Define the objective function
    function f_K_eqs(f, ξ, _)
        for i in 1 : rxn_system.n_reaction
            f_i1 = K_eqs[i]
            f_i2 = 1.0
            for j in 1 : rxn_system.n_species
                # Reactants
                if rxn_system.stoich[j, i] < 0
                    f_i1 *= (rxn_system.concs_init[j] + rxn_system.stoich[j, :] · ξ) ^ -rxn_system.stoich[j, i]
                # Products
                elseif rxn_system.stoich[j, i] > 0
                    f_i2 *= (rxn_system.concs_init[j] + rxn_system.stoich[j, :] · ξ) ^ rxn_system.stoich[j, i]
                end
            end
            f[i] = f_i1 - f_i2
        end
        nothing
    end
    # Compute reaction extents
    ξ = rxn_system.stoich \ (rxn_system.concs - rxn_system.concs_init)
    # Create objective function and Jacobian (automatic, with autodiff)
    problem = NonlinearSolve.NonlinearProblem(f_K_eqs, ξ; maxiters=maxiters, abstol=abstol, reltol=reltol)
    # Run the nonlinear optimization
    solution = NonlinearSolve.solve(problem, NonlinearSolve.TrustRegion())
    # If nonlinear optimization was successful, update `rxn_system.concs`
    if SciMLBase.successful_retcode(solution)
        rxn_system.concs .= rxn_system.concs_init
        rxn_system.concs .+= rxn_system.stoich * solution.u
    end
    # Return the nonlinear solution
    solution
end

function solve(rxn_system; φ=1.0, rate_consts=:forward, maxiters=1000, abstol=1.0e-9, reltol=0.0)
    if rate_consts === :forward
        # Evaluate K_eqs from forward rate constants and φ_fwd
        K_eqs = (rxn_system.fwd_rate_consts .- 2 * φ) ./ (rxn_system.fwd_rate_consts - φ)
    elseif rate_consts === :reverse
        # Evaluate K_eqs from reverse rate constants and φ_rev
        K_eqs = (rxn_system.rev_rate_consts .+ φ) ./ rxn_system.rev_rate_consts
    else
        throw(ArgumentError("invalid rate_consts value $rate_consts (must be :forward or :reverse)"))
    end
    solve(rxn_system, K_eqs; maxiters=maxiters, abstol=abstol, reltol=reltol)
end

end
