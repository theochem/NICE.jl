using NICE
using Test

@testset "NICE.jl" begin

    # Example 1

    test_concs = [0.9261879203, 0.9623865752, 0.0926187920]

    stoich = [[-0.5, 1.0, 0.0] [-0.5, -1.0, 1.0]]
    keq_vals = [1.0, 0.1]

    concs = [1.0, 0.2, 0.4]
    rxn_system = ReactionSystem(stoich, concs, keq_vals, φ=1.0)

    run_kmc(rxn_system, 1000000, ε=1.0e-5)
    @test isapprox(rxn_system.concs, test_concs, atol=1.0e-5)

    concs = [1.0, 0.2, 0.4]
    rxn_system = ReactionSystem(stoich, concs, keq_vals, φ=1.0)

    run_nekmc(rxn_system, 1000000, ε=1.0e-5)
    @test isapprox(rxn_system.concs, test_concs, atol=1.0e-5)

    # Example 2

    test_concs = [9.0e-01, 1.10212630e-08, 1.00002433e-10, 9.98999891e-02, 0.0, 9.99999e-05]

    stoich = [[-1, -1,  0,  1,  0,  0] [ 0, -1, -1,  0,  1,  0] [ 0,  0, -1, -1,  0,  1] [-1,  0,  0,  0, -1,  1]]
    keq_vals = [1.0e7, 1.0e9, 1.0e7, 1.0e9]

    concs = [1.0, 0.1, 1e-4, 0.0, 0.0, 0.0]
    rxn_system = ReactionSystem(stoich, concs, keq_vals, φ=1.0)

    run_kmc(rxn_system, 5000000, ε=1.0e-5)
    @test isapprox(rxn_system.concs, test_concs, atol=1.0e-5)

    concs = [1.0, 0.1, 1e-4, 0.0, 0.0, 0.0]
    rxn_system = ReactionSystem(stoich, concs, keq_vals, φ=1.0)

    run_nekmc(rxn_system, 1000000, ε=1.0e-5)
    @test isapprox(rxn_system.concs, test_concs, atol=1.0e-5)

end
