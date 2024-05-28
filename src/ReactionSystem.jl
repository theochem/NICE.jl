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
\\sum_{\\left. n \\middle| c_{1n} < 0 \\right.} &\\left|c_{1n}\\right| I_n
    &\\xrightleftharpoons[k_r^{\\left(1\\right)}]{k_f^{\\left(1\\right)}}
    \\sum_{\\left. n \\middle| c_{1n} > 0 \\right.} &\\left|c_{1n}\\right| I_n
    \\\\
\\vdots & ~ & ~ & \\vdots
    \\\\
\\sum_{\\left. n \\middle| c_{mn} < 0 \\right.} &\\left|c_{mn}\\right| I_n
    &\\xrightleftharpoons[k_r^{\\left(M\\right)}]{k_f^{\\left(M\\right)}}
    \\sum_{\\left. n \\middle| c_{Mn} > 0 \\right.} &\\left|c_{Mn}\\right| I_n
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
function ReactionSystem(
    stoich::AbstractMatrix{Float64},
    concs::AbstractVector{Float64},
    rev_rate_consts::AbstractVector{Float64},
    fwd_rate_consts::AbstractVector{Float64},
)
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
        n_species, n_reaction, stoich, rev_rate_consts, fwd_rate_consts,
        concs_init, n_iter, time, concs, rev_rates, fwd_rates, net_rates,
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
function ReactionSystem(
    stoich::AbstractMatrix{Float64},
    concs::AbstractVector{Float64},
    K_eqs::AbstractVector{Float64};
    φ::Real=1.0,
)
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
        n_species, n_reaction, stoich, rev_rate_consts, fwd_rate_consts,
        concs_init, n_iter, time, concs, rev_rates, fwd_rates, net_rates,
    )
end
