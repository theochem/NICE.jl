function update_rates(rxn_system::ReactionSystem, rev_stoich::Matrix{Float64}, fwd_stoich::Matrix{Float64})
  # r_i = k_i \prod_{j=1}^N c_j^{\max{0, +S_{ij}}}
  rxn_system.rev_rates .= rxn_system.rev_rate_consts
  for i = 1:rxn_system.n_reaction
    for j = 1:rxn_system.n_species
      rxn_system.rev_rates[i] *= rxn_system.concs[j] ^ rev_stoich[j, i]
    end
  end
  # f_i = k_i \prod_{j=1}^N c_j^{\max{0, -S_{ij}}}
  rxn_system.fwd_rates .= rxn_system.fwd_rate_consts
  for i = 1:rxn_system.n_reaction
    for j = 1:rxn_system.n_species
      rxn_system.fwd_rates[i] *= rxn_system.concs[j] ^ fwd_stoich[j, i]
    end
  end
  # n_i = f_i - r_i
  rxn_system.net_rates .= rxn_system.fwd_rates .- rxn_system.rev_rates
end

function initial_moving_average_rates(rxn_system::ReactionSystem, n_avg::Integer)
  # t = \sum_{i = 1}^M |n_i|
  t = 0
  for i = 1:rxn_system.n_reaction
    t += abs(rxn_system.net_rates[i])
  end
  return fill(t, n_avg)
end

function update_moving_average_rates_kmc(rxn_system::ReactionSystem, pvec::Vector{Float64}, tot_rates::Vector{Float64})
  # Rotate tot_rates, removing oldest entry
  circshift!(tot_rates, 1)
  # Insert new entry
  for i = 1:rxn_system.n_reaction
    tot_rates[1] += abs(rxn_system.net_rates[i])
  end
end

function update_moving_average_rates_nekmc(rxn_system::ReactionSystem, pvec::Vector{Float64}, tot_rates::Vector{Float64})
  # Rotate tot_rates, removing oldest entry
  circshift!(tot_rates, 1)
  # Insert new entry
  tot_rates[1] = pvec[end]
end

function update_step_size(
  rxn_system::ReactionSystem,
  tot_rates::Vector{Float64},
  tot_rate_prev::Float64,
  n_avg::Int,
  ε::Float64,
  ε_mult::Float64,
)
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
  return tot_rate_prev, ε
end

function update_pvec_kmc(rxn_system::ReactionSystem, pvec::Vector{Float64})
  # p_0 = 0
  # p_i = p_{i - 1} + |r_i|
  t = 0.
  for i = 1:rxn_system.n_reaction
    pvec[i] = t + abs(rxn_system.rev_rates[i])
    t = pvec[i]
  end
  # p_{M + i} = p_{M + i - 1} + |f_i|
  for i = 1:rxn_system.n_reaction
    pvec[rxn_system.n_reaction + i] = t + abs(rxn_system.fwd_rates[i])
    t = pvec[rxn_system.n_reaction + i]
  end
end

function update_pvec_nekmc(rxn_system::ReactionSystem, pvec::Vector{Float64})
  # p_0 = 0
  # p_i = p_{i - 1} + |n_i|
  t = 0.
  for i = 1:rxn_system.n_reaction
    pvec[i] = t + abs(rxn_system.net_rates[i])
    t = pvec[i]
  end
end

function do_random_reaction_kmc(
  rxn_system::ReactionSystem,
  pvec::Vector{Float64},
  concs::Vector{Float64},
  ε::Float64,
  ε0::Float64,
)
  # Select a reaction randomly
  @label select_reaction
  i_rxn = searchsortedfirst(pvec, pvec[end] * rand(Float64))

  # Do the reaction if it is possible
  concs .= rxn_system.concs
  if i_rxn <= rxn_system.n_reaction
    # Reverse reaction
    if rxn_system.rev_rates[i_rxn] > 0
      # c_{n + 1; j} = c_{n; j} - \sum_{k = 1}^N  \sgn(r_i) ε S_{ij}
      for i = 1:rxn_system.n_species
        concs[i] -= ε * rxn_system.stoich[i, i_rxn]
      end
    end
  else
    # Forward reaction
    # c_{n + 1; j} = c_{n; j} + \sum_{k = 1}^N \sgn(f_i) ε S_{ij}
    i_rxn -= rxn_system.n_reaction
    if rxn_system.fwd_rates[i_rxn] > 0
      for i = 1:rxn_system.n_species
        concs[i] += ε * rxn_system.stoich[i, i_rxn]
      end
    end
  end

  # Check that concentrations are >= 0
  clamp!(concs, 0, Inf)
  # for i = 1:rxn_system.n_species
  #   if concs[i] < 0
  #     # We need to re-try this, go back to the top
  #     @goto select_reaction
  #   end
  # end

  # Update the reaction system's concentrations
  rxn_system.concs .= concs

  # Update the time step
  # t = t - (ε / ε_0) \log(\rand()) / |p[2M]|
  rxn_system.time -= ε * log(rand(Float64)) / (ε0 * abs(pvec[end]))
end

function do_random_reaction_nekmc(
  rxn_system::ReactionSystem,
  pvec::Vector{Float64},
  concs::Vector{Float64},
  ε::Float64,
  ε0::Float64,
)
  # Select a reaction randomly
  @label select_reaction
  i_rxn = searchsortedfirst(pvec, pvec[end] * rand(Float64))

  # Do the reaction if it is possible
  # c_{n + 1; j} = c_{n; j} + \sum_{k = 1}^N \sgn(n_i) ε S_{ij}
  concs .= rxn_system.concs
  for i = 1:rxn_system.n_species
    concs[i] += ε * sign(rxn_system.net_rates[i_rxn]) * rxn_system.stoich[i, i_rxn]
  end
  # Check that concentrations are >= 0
  clamp!(concs, 0, Inf)
  # for i = 1:rxn_system.n_species
  #   if concs[i] < 0
  #     # We need to re-try this, go back to the top
  #     @goto select_reaction
  #   end
  # end

  # Update the reaction system's concentrations
  rxn_system.concs .= concs

  # Update the time step
  # t = t - (ε / ε_0) \log(\rand()) / |p[2M]|
  rxn_system.time -= ε * log(rand(Float64)) / (ε0 * abs(pvec[end]))
end

function collect_data_noop(rxn_system::ReactionSystem, i_data::Integer)
  return
end

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
  ε::Real=Float64(1.0e-3),
  ε_tol::Real=Float64(1.0e-12),
  ε_mult::Real=Float64(0.1),
  n_check::Integer=Int(100),
  n_avg::Integer=Int(100),
  n_data::Integer=Int(0),
  collect_data=collect_data_noop,
)
  # Working arrays
  rev_stoich = clamp.(+rxn_system.stoich, 0, Inf)
  fwd_stoich = clamp.(-rxn_system.stoich, 0, Inf)
  concs = zeros(Float64, rxn_system.n_species)
  pvec = zeros(Float64, 2 * rxn_system.n_reaction)

  # Update rates of reaction
  update_rates(rxn_system, rev_stoich, fwd_stoich)

  # Make vector of sum of abs(net_rates) for our moving average
  tot_rates = initial_moving_average_rates(rxn_system, n_avg)

  # Set the previous total rate
  tot_rate_prev = tot_rates[1]

  # Store original ε
  ε0 = ε

  # Collect initial simulation data point
  i_data = 1
  if n_data == 0 || collect_data === collect_data_noop
    n_data = typemax(Int)
  else
    collect_data(rxn_system, i_data)
    i_data += 1
  end

  # Begin iterating
  rxn_system.n_iter = 0
  # Stop if number of iterations exceeded or ε converged
  while (rxn_system.n_iter != n_iter) && (ε > ε_tol)
    # Compute cumulative sum of abs(net_rates)
    update_pvec_kmc(rxn_system, pvec)

    # Update total rate moving average
    update_moving_average_rates_kmc(rxn_system, pvec, tot_rates)

    # Check every n_check iterations for a significant (in|de)crease in total rate moving average
    if (rxn_system.n_iter % n_check == 0)
      tot_rate_prev, ε = update_step_size(rxn_system, tot_rates, tot_rate_prev, n_avg, ε, ε_mult)
    end

    # Do random reaction
    do_random_reaction_kmc(rxn_system, pvec, concs, ε, ε0)

    # Update rates of reaction
    update_rates(rxn_system, rev_stoich, fwd_stoich)

    # Collect simulation data
    if rxn_system.n_iter % n_data == 0
      collect_data(rxn_system, i_data)
      i_data += 1
    end

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
  ε::Real=Float64(1.0e-3),
  ε_tol::Real=Float64(1.0e-12),
  ε_mult::Real=Float64(0.1),
  n_check::Integer=Int(100),
  n_avg::Integer=Int(100),
  n_data::Integer=Int(0),
  collect_data=collect_data_noop,
)
  # Working arrays
  rev_stoich = clamp.(+rxn_system.stoich, 0, Inf)
  fwd_stoich = clamp.(-rxn_system.stoich, 0, Inf)
  concs = zeros(Float64, rxn_system.n_species)
  pvec = zeros(Float64, rxn_system.n_reaction)

  # Update rates of reaction
  update_rates(rxn_system, rev_stoich, fwd_stoich)

  # Make vector of sum of abs(net_rates) for our moving average
  tot_rates = initial_moving_average_rates(rxn_system, n_avg)

  # Set the previous total rate
  tot_rate_prev = tot_rates[1]

  # Store original ε
  ε0 = ε

  # Collect initial simulation data point
  i_data = 1
  if n_data == 0 || collect_data === collect_data_noop
    n_data = typemax(Int)
  else
    collect_data(rxn_system, i_data)
    i_data += 1
  end

  # Begin iterating
  rxn_system.n_iter = 0
  # Stop if number of iterations exceeded or ε converged
  while (rxn_system.n_iter != n_iter) && (ε > ε_tol)
    # Compute cumulative sum of abs(net_rates)
    update_pvec_nekmc(rxn_system, pvec)

    # Update total rate moving average
    update_moving_average_rates_nekmc(rxn_system, pvec, tot_rates)

    # Check every n_check iterations for a significant (in|de)crease in total rate moving average
    if (rxn_system.n_iter % n_check == 0)
      tot_rate_prev, ε = update_step_size(rxn_system, tot_rates, tot_rate_prev, n_avg, ε, ε_mult)
    end

    # Do random reaction
    do_random_reaction_nekmc(rxn_system, pvec, concs, ε, ε0)

    # Update rates of reaction
    update_rates(rxn_system, rev_stoich, fwd_stoich)

    # Collect simulation data point
    if rxn_system.n_iter % n_data == 0
      collect_data(rxn_system, i_data)
      i_data += 1
    end

    # Go to next iteration
    rxn_system.n_iter += 1
  end
end
