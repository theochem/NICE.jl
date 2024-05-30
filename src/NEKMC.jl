using LinearAlgebra

export simulate!

"""
    simulate!(rxn_system;
             n_iter=Int(1e+8), chunk_iter=Int(1e+4),
             ε=1.0e-4, ε_scale=1.0, ε_concs=0.0,
             tol_t=Inf, tol_ε=0.0, tol_concs=0.0)

Run a *N*et-*E*vent *K*inetic *M*onte *C*arlo (NEKMC) simulation to find the
equilibrium concentrations of the reaction system.
"""
function simulate!(
    rxn_system::ReactionSystem;
    n_iter::Int=Int(1e+8),
    chunk_iter::Int=Int(1e+4),
    ε::Float64=1.0e-4,
    ε_scale::Float64=1.0,
    ε_concs::Float64=0.0,
    tol_t::Float64=Inf,
    tol_ε::Float64=0.0,
    tol_concs::Float64=0.0,
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

function update_rates(
    rxn_system::ReactionSystem,
)
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

function select_reaction(
    rxn_system::ReactionSystem,
    pvec::Vector{Float64},
)
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

function do_reaction(
    rxn_system::ReactionSystem,
    i_rxn::Int,
    ε::Float64,
)
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
