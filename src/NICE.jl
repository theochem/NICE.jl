module NICE


export ReactionSystem, run_kmc, run_nekmc


mutable struct ReactionSystem

    # Constants
    n_reaction::Int
    n_species::Int
    stoich::Matrix{Float64}
    rev_rate_consts::Vector{Float64}
    fwd_rate_consts::Vector{Float64}

    # Mutables
    n_iter::Int
    time::Float64
    concs::Vector{Float64}
    rev_rates::Vector{Float64}
    fwd_rates::Vector{Float64}
    net_rates::Vector{Float64}
end


function ReactionSystem(stoich, concs, rev_rate_consts, fwd_rate_consts)
    n_reaction = size(stoich, 1)
    n_species = size(stoich, 2)
    n_iter = 0
    time = 0.
    rev_rates = zeros(Float64, n_reaction)
    fwd_rates = zeros(Float64, n_reaction)
    net_rates = zeros(Float64, n_reaction)
    return ReactionSystem(n_reaction, n_species, stoich, rev_rate_consts, fwd_rate_consts,
        n_iter, time, concs, rev_rates, fwd_rates, net_rates)
end


function ReactionSystem(stoich, concs, keq_vals, φ::Real=1.0)
    n_reaction = size(stoich, 1)
    n_species = size(stoich, 2)
    rev_rate_consts = φ ./ (keq_vals .+ 1)
    fwd_rate_consts = φ .- rev_rate_consts
    n_iter = 0
    time = 0.
    rev_rates = zeros(Float64, n_reaction)
    fwd_rates = zeros(Float64, n_reaction)
    net_rates = zeros(Float64, n_reaction)
    return ReactionSystem(n_reaction, n_species, stoich, rev_rate_consts, fwd_rate_consts,
        n_iter, time, concs, rev_rates, fwd_rates, net_rates)
end


function run_kmc(rxn_system, n_iter, ε=1.0e-4)
    pvec = zeros(2 * rxn_system.n_reaction)
    for i in 1 : n_iter
        update_rates(rxn_system)
        i_rxn = kmc_select_reaction(rxn_system, pvec)
        rxn_system.time += kmc_do_reaction(rxn_system, i_rxn, ε)
        rxn_system.n_iter += 1
    end
end


function run_nekmc(rxn_system, n_iter, ε=1.0e-4)
    pvec = zeros(rxn_system.n_reaction)
    for i in 1 : n_iter
        update_rates(rxn_system)
        i_rxn = nemkc_select_reaction(rxn_system, pvec)
        rxn_system.time += nemkc_do_reaction(rxn_system, i_rxn, ε)
        rxn_system.n_iter += 1
    end
end


function update_rates(rxn_system)
    for i in 1 : rxn_system.n_reaction
        rev_rate = rxn_system.rev_rate_consts[i]
        fwd_rate = rxn_system.fwd_rate_consts[i]
        for j in 1 : rxn_system.n_species
            s = rxn_system.stoich[i, j]
            if s >= 0.
                rev_rate *= rxn_system.concs[j] ^ s
            else
                fwd_rate *= rxn_system.concs[j] ^ abs(s)
            end
        end
        rxn_system.rev_rates[i] = rev_rate
        rxn_system.fwd_rates[i] = fwd_rate
        rxn_system.net_rates[i] = fwd_rate - rev_rate
    end
end


function kmc_select_reaction(rxn_system, pvec)
    # Update probability vector
    p = 0.
    for i in 1 : rxn_system.n_reaction
        p += rxn_system.fwd_rates[i]
        pvec[i] = p
    end
    for i in 1 : rxn_system.n_reaction
        p += rxn_system.rev_rates[i]
        pvec[i + rxn_system.n_reaction] = p
    end
    # Select random reaction
    p *= rand(Float64)
    for i in 1 : 2 * rxn_system.n_reaction
        if pvec[i] > p
            return i
        end
    end
    # sentinel return value
    return 2 * rxn_system.n_reaction
end


function nekmc_select_reaction(rxn_system, pvec)
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
    return rxn_system.n_reaction
end


function kmc_do_reaction(rxn_system, i_rxn, ε)
    # Forward reaction (i_rxn <= n_rxn)
    if i_rxn <= rxn_system.n_reaction
        rate = rxn_system.fwd_rates[i_rxn]
        if rate >= 0
            for j = 1 : rxn_system.n_species
                rxn_system.concs[j] += rxn_system.stoich[i_rxn, j] * ε
            end
            return ε / rate
        end
    # Reverse reaction (i_rxn > n_rxn)
    else
        i_rxn -= rxn_system.n_reaction
        rate = rxn_system.rev_rates[i_rxn]
        if rate >= 0
            for j = 1 : rxn_system.n_species
                rxn_system.concs[j] -= rxn_system.stoich[i_rxn, j] * ε
            end
            return ε / rate
        end
    end
    # Sentinel return value
    return 0.
end


function nekmc_do_reaction(rxn_system, i_rxn, ε)
    # Net reaction
    rate = rxn_system.net_rates[i_rxn]
    if rate >= 0
        for j = 1 : rxn_system.n_species
            rxn_system.concs[j] += rxn_system.stoich[i_rxn, j] * ε
        end
        return ε / rate
    else
        for j = 1 : rxn_system.n_species
            rxn_system.concs[j] -= rxn_system.stoich[i_rxn, j] * ε
        end
        return -ε / rate
    end
end


end
