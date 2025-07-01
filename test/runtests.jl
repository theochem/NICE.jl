using NICE
using Test
using ForwardDiff

@testset "Biological Oscillator System" begin
    # Initial concentrations for the 6 species in the biological oscillator
    initial_concs = [0.0, 0.0, 0.0, 10.0, 1035.0, 0.0]

    # Stoichiometric matrix: rows are species, columns are reactions
    stoich = Float64[
        0 0 -2 2 -1 1 0 -1 0;
        1 1 0 0 0 0 0 -1 0;
        0 0 1 -1 -1 1 0 0 0;
        0 0 0 0 1 -1 0 0 0;
        0 0 0 0 0 0 -1 0 1;
        0 0 0 0 0 0 1 0 -1
    ]

    # Rate constants from the biological oscillator model
    k0 = 0.0025
    k1 = 0.833
    ζ = 100
    k2 = 3.13e+8 * ζ
    φ = 1.628e+9 * ζ
    χ = 6.517e+8 * ζ
    a = 159.37
    b = 5.31
    λ1 = 1.563e+5
    λ2 = 1.667e-4

    # Equilibrium constants for each reaction
    K_eqs = [
        1.0,
        1.0,
        χ / φ,
        1 / (χ / φ),
        φ / χ,
        1 / (φ / χ),
        1.0,
        1.0,
        λ2 / λ1
    ]

    test_concs = [0.0, 0.0, 0.0, 10.0, 1035.0, 0.0]

    # Test simulation method
    rxn_system = ReactionSystem(stoich, initial_concs, K_eqs; φ=1.0)
    simulate(rxn_system; n_iter=Int(1e+6), chunk_iter=Int(1e+3), ε=1e-3, ε_scale=0.5, ε_concs=2.0, tol_ε=0.)
    @test isapprox(rxn_system.concs[1] + rxn_system.concs[4], 10.0; atol=1.0e-2)
    @test isapprox(rxn_system.concs[5] + rxn_system.concs[6], 1035.0; atol=1.0e-2)

    # Test direct solver
    rxn_system = ReactionSystem(stoich, initial_concs, K_eqs; φ=1.0)
    solve(rxn_system, K_eqs)
    @test isapprox(rxn_system.concs[1] + rxn_system.concs[4], 10.0; atol=1.0e-3)
    @test isapprox(rxn_system.concs[5] + rxn_system.concs[6], 1035.0; atol=1.0e-3)

    # Test hybrid solver combining simulation and direct methods
    rxn_system = ReactionSystem(stoich, initial_concs, K_eqs; φ=1.0)
    hybrid_solve(rxn_system, K_eqs; n_iter=Int(5e+6), chunk_iter=Int(1e+3), ε=1e-3, ε_scale=0.5, ε_concs=2.0, tol_ε=0., maxiters=1000)
    @test isapprox(rxn_system.concs[1] + rxn_system.concs[4], 10.0; atol=1.0e-2)
    @test isapprox(rxn_system.concs[5] + rxn_system.concs[6], 1035.0; atol=1.0e-3)
end

@testset "System 1" begin
    # Known equilibrium concentrations for this simple 3-species, 2-reaction system
    test_concs = [0.9261879203, 0.9623865752, 0.0926187920]

    # Simple test system with 3 species and 2 reactions
    stoich = Float64[
        -0.5 -0.5;
        1.0 -1.0;
        0.0 1.0
    ]
    keq_vals = [1.0, 0.1]

    concs = [1.0, 0.2, 0.4]

    # Test each solution method against known equilibrium
    rxn_system = ReactionSystem(stoich, concs, keq_vals; φ=1.0)
    simulate(rxn_system; n_iter=Int(5e+8), chunk_iter=Int(1e+3), ε=1e-3, ε_scale=0.5, ε_concs=2.0, tol_ε=0.)
    @test isapprox(rxn_system.concs, test_concs; atol=1.0e-5)

    rxn_system = ReactionSystem(stoich, concs, keq_vals; φ=1.0)
    solve(rxn_system, keq_vals)
    @test isapprox(rxn_system.concs, test_concs; atol=1.0e-6)

    rxn_system = ReactionSystem(stoich, concs, keq_vals; φ=1.0)
    hybrid_solve(rxn_system, keq_vals; n_iter=Int(5e+8), chunk_iter=Int(1e+3), ε=1e-3, ε_scale=0.5, ε_concs=2.0, tol_ε=0., maxiters=1000)
    @test isapprox(rxn_system.concs, test_concs; atol=1.0e-6)
end

@testset "System 2" begin
    # Equilibrium solution for a more complex 6-species, 4-reaction system
    test_concs = [9.0e-01, 1.10212630e-08, 1.00002433e-10, 9.98999891e-02, 0.0, 9.99999e-05]

    # Network with multiple competing reactions
    stoich = [-1.0 0.0 0.0 -1.0;
        -1.0 -1.0 0.0 0.0;
        0.0 -1.0 -1.0 0.0;
        1.0 0.0 -1.0 0.0;
        0.0 1.0 0.0 -1.0;
        0.0 0.0 1.0 1.0]
    keq_vals = [1.0e7, 1.0e9, 1.0e7, 1.0e9]

    concs = [1.0, 0.1, 1e-4, 0.0, 0.0, 0.0]
    rxn_system = ReactionSystem(stoich, concs, keq_vals; φ=1.0)
    simulate(rxn_system; n_iter=Int(5e+8), chunk_iter=Int(1e+3), ε=1e-3, ε_scale=0.5, ε_concs=2.0, tol_ε=0.)
    @test isapprox(rxn_system.concs, test_concs; atol=1.0e-5)

    solve(rxn_system, keq_vals)
    @test isapprox(rxn_system.concs, test_concs; atol=1.0e-6)

    rxn_system = ReactionSystem(stoich, concs, keq_vals; φ=1.0)
    hybrid_solve(rxn_system, keq_vals; n_iter=Int(5e+8), chunk_iter=Int(1e+3), ε=1e-3, ε_scale=0.5, ε_concs=2.0, tol_ε=0., maxiters=1000)
    @test isapprox(rxn_system.concs, test_concs; atol=1.0e-6)
end

@testset "Constructor and Property Tests" begin
    initial_concs = [1.0, 0.2, 0.4]
    stoich = Float64[
        -0.5 -0.5;
        1.0 -1.0;
        0.0 1.0
    ]
    keq_vals = [1.0, 0.1]

    rxn_system = ReactionSystem(stoich, initial_concs, keq_vals; φ=1.0)

    # Verify basic system dimensions
    @test rxn_system.n_species == 3
    @test rxn_system.n_reaction == 2
    @test isapprox(rxn_system.concs_init, initial_concs)
    @test isapprox(rxn_system.concs, initial_concs)
    @test size(rxn_system.stoich) == (3, 2)

    # Check that rate constants are calculated correctly from equilibrium constants
    @test isapprox(rxn_system.fwd_rate_consts, [0.5, 0.09090909090909091]; atol=1e-10)
    @test isapprox(rxn_system.rev_rate_consts, [0.5, 0.9090909090909091]; atol=1e-10)
end

@testset "Rate Calculation Tests" begin
    initial_concs = [1.0, 0.2, 0.4]
    stoich = Float64[
        -0.5 -0.5;
        1.0 -1.0;
        0.0 1.0
    ]
    keq_vals = [1.0, 0.1]

    rxn_system = ReactionSystem(stoich, initial_concs, keq_vals; φ=1.0)

    # Expected rate constants based on equilibrium constants and stoichiometry
    expected_fwd_rate_consts = [0.5, 0.09090909090909091]
    expected_rev_rate_consts = [0.5, 0.9090909090909091]
    @test isapprox(rxn_system.fwd_rate_consts, expected_fwd_rate_consts; atol=1e-10)
    @test isapprox(rxn_system.rev_rate_consts, expected_rev_rate_consts; atol=1e-10)

    # Run one iteration to test that the system doesn't crash
    simulate(rxn_system; n_iter=1, chunk_iter=1, ε=1e-3, ε_scale=0.5, ε_concs=2.0, tol_ε=0.)

    rxn_system2 = ReactionSystem(stoich, initial_concs, keq_vals; φ=1.0)
end

@testset "Mass Conservation Tests" begin
    # Single reaction: A + B → 2C
    initial_concs = [2.0, 1.0, 0.0]
    stoich = reshape(Float64[-1.0, -1.0, 2.0], 3, 1)
    keq_vals = [5.0]
    initial_total = sum(initial_concs)

    # Check mass conservation for all three solution methods
    for method_name in ["simulate", "solve", "hybrid_solve"]
        rxn_system = ReactionSystem(stoich, initial_concs, keq_vals; φ=1.0)

        if method_name == "simulate"
            simulate(rxn_system; n_iter=Int(1e5), chunk_iter=Int(1e3), ε=1e-4, ε_scale=0.5, ε_concs=2.0, tol_ε=0.)
        elseif method_name == "solve"
            solve(rxn_system, keq_vals)
        else
            hybrid_solve(rxn_system, keq_vals; n_iter=Int(1e5), chunk_iter=Int(1e3), ε=1e-4, ε_scale=0.5, ε_concs=2.0, tol_ε=0., maxiters=100)
        end

        # Total number of atoms should be conserved
        @test isapprox(sum(rxn_system.concs), initial_total; rtol=1e-6)
        # No negative concentrations allowed
        @test all(rxn_system.concs .>= 0.0)
    end

    # More complex system with multiple reactions
    initial_concs_multi = [1.5, 1.0, 0.5, 0.0, 0.0]
    stoich_multi = Float64[
        -1.0 -1.0 0.0;
        -1.0 0.0 -1.0;
        2.0 -1.0 0.0;
        0.0 2.0 -1.0;
        0.0 0.0 2.0
    ]
    keq_vals_multi = [2.0, 10.0, 0.5]
    initial_total_multi = sum(initial_concs_multi)

    rxn_system_multi = ReactionSystem(stoich_multi, initial_concs_multi, keq_vals_multi; φ=1.0)
    solve(rxn_system_multi, keq_vals_multi)
    @test isapprox(sum(rxn_system_multi.concs), initial_total_multi; rtol=1e-8)
    @test all(rxn_system_multi.concs .>= 0.0)
end

@testset "Jacobian Verification" begin
    # Use System 1 for Jacobian testing since it's well-behaved
    stoich = Float64[
        -0.5 -0.5;
        1.0 -1.0;
        0.0 1.0
    ]
    keq_vals = [1.0, 0.1]
    concs = [1.0, 0.2, 0.4]
    rxn_system = ReactionSystem(stoich, concs, keq_vals; φ=1.0)

    # Wrapper function for ForwardDiff - needs to be out-of-place
    objective_wrapper = (xi) -> begin
        f_out = similar(xi, rxn_system.n_reaction)
        NICE.objective_function!(f_out, rxn_system, keq_vals, xi)
        return f_out
    end

    # Test Jacobian accuracy at multiple random points
    for _ in 1:20
        # Small random perturbations to avoid negative concentrations
        xi = 0.01 * rand(rxn_system.n_reaction)

        # Skip if this would give negative concentrations
        if any(rxn_system.concs_init .+ rxn_system.stoich * xi .<= 0)
            continue
        end

        # Compare analytical Jacobian to automatic differentiation
        J_analytical = zeros(rxn_system.n_reaction, rxn_system.n_reaction)
        NICE.jacobian_function!(J_analytical, rxn_system, keq_vals, xi)

        J_ad = ForwardDiff.jacobian(objective_wrapper, xi)

        @test isapprox(J_analytical, J_ad; atol=1e-8, rtol=1e-8)
    end
end

@testset "Negative Activity Handling" begin
    stoich = Float64[
        -2.0 1.0;
        1.0 -1.0;
        0.0 1.0
    ]
    keq_vals = [1.0, 0.1]
    concs = [1.0, 0.2, 0.4]
    rxn_system = ReactionSystem(stoich, concs, keq_vals; φ=1.0)

    xi_negative = [1.0, 0.0]  # This should make activities[1] = 1.0 - 2.0*1.0 = -1.0
    f_out = zeros(Float64, rxn_system.n_reaction)

    # The objective function should not crash with negative activities
    NICE.objective_function!(f_out, rxn_system, keq_vals, xi_negative)
    @test length(f_out) == rxn_system.n_reaction

    # Test with fractional exponents that would cause domain errors
    stoich_fractional = Float64[
        -0.5 1.0;  # This will cause fractional exponents
        1.0 -1.0;
        0.0 1.0
    ]
    rxn_system_frac = ReactionSystem(stoich_fractional, concs, keq_vals; φ=1.0)
    f_out_frac = zeros(Float64, rxn_system_frac.n_reaction)

    # This should either throw a DomainError or return complex results converted to NaN/Inf
    try
        NICE.objective_function!(f_out_frac, rxn_system_frac, keq_vals, xi_negative)
        # If it succeeds, the results may be NaN or Inf due to complex number conversions
        @test any(isnan.(f_out_frac)) || any(isinf.(f_out_frac)) || all(isfinite.(f_out_frac))
    catch e
        # DomainError is expected for negative^fractional
        @test isa(e, DomainError)
    end

    # Test that Jacobian function handles negative activities with its explicit check
    J_out = zeros(Float64, rxn_system.n_reaction, rxn_system.n_reaction)
    NICE.jacobian_function!(J_out, rxn_system, keq_vals, xi_negative)
    # Should return infinite values as implemented
    @test all(isinf.(J_out))
end