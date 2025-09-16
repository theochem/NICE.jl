"""
NICE (N-species Ice Table) module -- provides simultaneous equilibria solvers:
- [`kmc_simulate`](@ref) kinetic Monte Carlo solver
- [`nekmc_simulate`](@ref) Net-event kinetic Monte Carlo solver
- [`solve`](@ref) Nonlinear equations solver (Newton trust region method)

The [`ReactionSystem`](@ref) type is provided as an interface for these solvers.
"""
module NICE

export ReactionSystem
export kmc_simulate
export nekmc_simulate
export solve
export kmc_hybrid_solve
export nekmc_hybrid_solve

include("ReactionSystem.jl")
include("KMC.jl")
include("Exact.jl")

end
