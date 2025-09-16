using LinearAlgebra

"""
    kmc_simulate(rxn_system;
             n_iter=Int(1e+8), ε=1.0e-4, ε_tol=1.0e-12, ε_mult=0.1,
             n_check=100, n_avg=100)

Run a *K*inetic *M*onte *C*arlo (NEKMC) simulation to find the
equilibrium concentrations of the reaction system.
"""
function kmc_simulate(
  rxn_system::ReactionSystem;
  n_iter::Integer=Int(1e+8),
  ε::Real=1.0e-3,
  ε_tol::Real=1.0e-12,
  ε_mult::Real=0.1,
  n_check=100,
  n_avg=100,
)
  # Working arrays
  rev_stoich = clamp.(+rxn_system.stoich, 0, Inf)
  fwd_stoich = clamp.(-rxn_system.stoich, 0, Inf)
  concs = zeros(Float64, rxn_system.n_species)

  # Update rates of reaction
  rxn_system.rev_rates .= rxn_system.rev_rate_consts .* vec(prod(rxn_system.concs .^ rev_stoich; dims=1))
  rxn_system.fwd_rates .= rxn_system.fwd_rate_consts .* vec(prod(rxn_system.concs .^ fwd_stoich; dims=1))
  rxn_system.net_rates .= rxn_system.fwd_rates .- rxn_system.rev_rates

  # Make vector of sum of abs(net_rates) for our moving average
  tot_rates = fill(sum(abs.(rxn_system.rev_rates) + abs.(rxn_system.fwd_rates)), n_avg)

  # Set the previous total rate
  tot_rate_prev = tot_rates[1]

  # Store original ε
  ε0 = ε

  # Begin iterating
  rxn_system.n_iter = 0
  # Stop if number of iterations exceeded or ε converged
  while (rxn_system.n_iter != n_iter) && (ε > ε_tol)

    # Compute cumulative sum of abs(net_rates)
    pvec = cumsum(cat(abs.(rxn_system.rev_rates), abs.(rxn_system.fwd_rates); dims=1); dims=1)

    # Update total rate moving average
    circshift!(tot_rates, 1)
    tot_rates[1] = pvec[end]

    # Check every n_check iterations for a significant (in|de)crease in total rate moving average
    if (n_iter % n_check == 0)
      # Compute total rate from moving average
      tot_rate_curr = sum(tot_rates) / n_avg
      if (tot_rate_curr < tot_rate_prev * ε_mult)
        # Set new previous total rate
        tot_rate_prev = tot_rate_curr

        # Scale ε
        ε *= ε_mult
      elseif (tot_rate_curr * ε_mult > tot_rate_prev)
        # Set new previous total rate
        tot_rate_prev = tot_rate_curr

        # Scale ε
        ε /= ε_mult
      end
    end

    # Select a reaction randomly
    i_rxn = searchsortedfirst(pvec, pvec[end] * rand(Float64))

    # Do the reaction if it is possible
    if i_rxn <= rxn_system.n_reaction
      concs .= rxn_system.concs .- (ε * sign(rxn_system.rev_rates[i_rxn])) .* rxn_system.stoich[:, i_rxn]
    else
      i_rxn -= rxn_system.n_reaction
      concs .= rxn_system.concs .+ (ε * sign(rxn_system.fwd_rates[i_rxn])) .* rxn_system.stoich[:, i_rxn]
    end
    any(concs .< 0) && continue
    rxn_system.concs .= concs

    # Update the time step
    rxn_system.time -= ε * log(rand(Float64)) / (ε0 * abs(pvec[end]))

    # Update rates of reaction
    rxn_system.rev_rates .= rxn_system.rev_rate_consts .* vec(prod(rxn_system.concs .^ rev_stoich; dims=1))
    rxn_system.fwd_rates .= rxn_system.fwd_rate_consts .* vec(prod(rxn_system.concs .^ fwd_stoich; dims=1))
    rxn_system.net_rates .= rxn_system.fwd_rates .- rxn_system.rev_rates

    # Go to next iteration
    rxn_system.n_iter += 1
  end
end

"""
    nekmc_simulate(rxn_system;
             n_iter=Int(1e+8), ε=1.0e-4, ε_tol=1.0e-12, ε_mult=0.1,
             n_check=100, n_avg=100)

Run a *N*et-*E*vent *K*inetic *M*onte *C*arlo (NEKMC) simulation to find the
equilibrium concentrations of the reaction system.
"""
function nekmc_simulate(
  rxn_system::ReactionSystem;
  n_iter::Integer=Int(1e+8),
  ε::Real=1.0e-3,
  ε_tol::Real=1.0e-12,
  ε_mult::Real=0.1,
  n_check=100,
  n_avg=100,
)
  # Working arrays
  rev_stoich = clamp.(+rxn_system.stoich, 0, Inf)
  fwd_stoich = clamp.(-rxn_system.stoich, 0, Inf)
  concs = zeros(Float64, rxn_system.n_species)
  pvec = zeros(Float64, rxn_system.n_reaction)

  # Update rates of reaction
  rxn_system.rev_rates .= rxn_system.rev_rate_consts .* vec(prod(rxn_system.concs .^ rev_stoich; dims=1))
  rxn_system.fwd_rates .= rxn_system.fwd_rate_consts .* vec(prod(rxn_system.concs .^ fwd_stoich; dims=1))
  rxn_system.net_rates .= rxn_system.fwd_rates .- rxn_system.rev_rates

  # Make vector of sum of abs(net_rates) for our moving average
  tot_rates = fill(sum(abs.(rxn_system.net_rates)), n_avg)

  # Set the previous total rate
  tot_rate_prev = tot_rates[1]

  # Store original ε
  ε0 = ε

  # Begin iterating
  rxn_system.n_iter = 0
  # Stop if number of iterations exceeded or ε converged
  while (rxn_system.n_iter != n_iter) && (ε > ε_tol)

    # Compute cumulative sum of abs(net_rates)
    cumsum!(pvec, abs.(rxn_system.net_rates))

    # Update total rate moving average
    circshift!(tot_rates, 1)
    tot_rates[1] = pvec[end]

    # Check every n_check iterations for a significant (in|de)crease in total rate moving average
    if (n_iter % n_check == 0)
      # Compute total rate from moving average
      tot_rate_curr = sum(tot_rates) / n_avg
      if (tot_rate_curr < tot_rate_prev * ε_mult)
        # Set new previous total rate
        tot_rate_prev = tot_rate_curr

        # Scale ε
        ε *= ε_mult
      elseif (tot_rate_curr * ε_mult > tot_rate_prev)
        # Set new previous total rate
        tot_rate_prev = tot_rate_curr

        # Scale ε
        ε /= ε_mult
      end
    end

    # Select a reaction randomly
    i_rxn = searchsortedfirst(pvec, pvec[end] * rand(Float64))

    # Do the reaction if it is possible
    concs .= rxn_system.concs .+ (ε * sign(rxn_system.net_rates[i_rxn])) .* rxn_system.stoich[:, i_rxn]
    any(concs .< 0) && continue
    rxn_system.concs .= concs

    # Update the time step
    rxn_system.time -= ε * log(rand(Float64)) / (ε0 * abs(pvec[end]))

    # Update rates of reaction
    rxn_system.rev_rates .= rxn_system.rev_rate_consts .* vec(prod(rxn_system.concs .^ rev_stoich; dims=1))
    rxn_system.fwd_rates .= rxn_system.fwd_rate_consts .* vec(prod(rxn_system.concs .^ fwd_stoich; dims=1))
    rxn_system.net_rates .= rxn_system.fwd_rates .- rxn_system.rev_rates

    # Go to next iteration
    rxn_system.n_iter += 1
  end
end
