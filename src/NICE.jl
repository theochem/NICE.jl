module NICE


export ReactionSystem
# export run_kmc, run_nekmc


# Define the struct with information pertaining to the reaction system
struct ReactionSystem
    n_reaction::Int
    n_species::Int
    stoich::Matrix{Float64}
    rev_rate_consts::Vector{Float64}
    fwd_rate_consts::Vector{Float64}
end


# Constructor #1
function ReactionSystem(stoich, rev_rate_consts, fwd_rate_consts::AbstractVector{Real})
    n_reaction = size(stoich, 1)
    n_species = size(stoich, 2)
    return ReactionSystem(n_reaction, n_species, stoich, rev_rate_consts, fwd_rate_consts)
end


# Constructor #2
function ReactionSystem(stoich, keq_vals, φ::Real=1.0)
    n_reaction = size(stoich, 1)
    n_species = size(stoich, 2)
    rev_rate_consts = φ ./ (keq_vals .+ 1)
    fwd_rate_consts = φ .- rev_rate_consts
    return ReactionSystem(n_reaction, n_species, stoich, rev_rate_consts, fwd_rate_consts)
end


# Functions for KMC and NEKMC. There are two definitions for each function name;
# one for KMC and one for NEKMC. They have different arguments!


function update_rates(rxn_system, concs, rev_rates, fwd_rates, net_rates)
    for i in 1 : rxn_system.n_reaction
        p = rxn_system.rev_rate_consts[i]
        q = rxn_system.fwd_rate_consts[i]
        for j in 1 : rxn_system.n_species
            s = rxn_system.stoich[i, j]
            if s >= 0.
                p *= concs[j] ^ s
            else
                q *= concs[j] ^ abs(s)
            end
        end
        rev_rates[i] = p
        fwd_rates[i] = q
        net_rates[i] = q - p
    end
end


function select_reaction(rxn_system, rev_rates, fwd_rates)
    #
    # ...
    #
    return 1 # index of reaction in ``stoich`` matrix.
end


function select_reaction(rxn_system, net_rates)
    #
    # ...
    #
    return 1 # index of reaction in ``stoich`` matrix.
end


function do_reaction(rxn_system, concs, rev_rates, fwd_rates, ε)
    #
    # Update the ``concs`` vector.
    # ε is the concentration step size.
    #
end


function do_reaction(rxn_system, concs, net_rates, ε)
    #
    # Update the ``concs`` vector.
    # ε is the concentration step size.
    #
end


# And now write the entire (NE)KMC algorithm from these pieces:
#
# function run_(ne)kmc( ..., ..., ..., ... )
#     ...
#     ...
#     ...
# end


end
