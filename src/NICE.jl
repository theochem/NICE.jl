"""
NICE (N-species Ice Table) module -- provides simultaneous equilibria solvers:
- [`simulate`](@ref) Net-event kinetic Monte Carlo solver
- [`solve`](@ref) Nonlinear equations solver (Newton trust region method)

The [`ReactionSystem`](@ref) type is provided as an interface for these solvers.
"""
module NICE

export ReactionSystem
export simulate
export solve

include("ReactionSystem.jl")
include("NEKMC.jl")
include("Exact.jl")

end # module NICE
