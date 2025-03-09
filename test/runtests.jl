using NICE
using Test

struct TestSystem
    stoich::Matrix{Float64}
    K_eqs::Vector{Float64}
    concs::Vector{Float64}
    concs_eq::Vector{Float64}
end

system1 = TestSystem(
    # stoich
    [-0.5 -0.5;
      1.0 -1.0;
      0.0  1.0],
    # K_eqs
    [1.0, 0.1],
    # concs
    [1.0, 0.2, 0.4],
    # concs_eq
    [0.9261879203, 0.9623865752, 0.0926187920],
)

system2 = TestSystem(
    # stoich
    [-1.0  0.0  0.0 -1.0;
     -1.0 -1.0  0.0  0.0;
      0.0 -1.0 -1.0  0.0;
      1.0  0.0 -1.0  0.0;
      0.0  1.0  0.0 -1.0;
      0.0  0.0  1.0  1.0],
    # K_eqs
    [1.0e7, 1.0e9, 1.0e7, 1.0e9],
    # concs
    [1.0, 0.1, 1e-4, 0.0, 0.0, 0.0],
    # concs_eq
    [9.0e-01, 1.10212630e-08, 1.00002433e-10, 9.98999891e-02, 0.0, 9.99999e-05],
)

@testset verbose=true "System 1" begin
    # Set up reaction system
    rxn_system1 = ReactionSystem(system1.stoich, copy(system1.concs), system1.K_eqs; φ=1.0)
    # Run and test NEKMC
    simulate!(
        rxn_system1;
        n_iter=Int(1e+9),
        chunk_iter=Int(1e+3),
        ε=1e-3,
        ε_scale=0.5,
        ε_concs=2.0,
        tol_ε=0.0,
    )
    @test isapprox(rxn_system1.concs, system1.concs_eq; rtol=1.0e-4)
    # Run and test exact nonlinear solver
    rxn_system1_copy = deepcopy(rxn_system1)
    solve!(rxn_system1_copy, system1.K_eqs; reltol=1.0e-6)
    @test isapprox(rxn_system1_copy.concs, system1.concs_eq; rtol=1.0e-6)
    # Run and test CMA solver
    rxn_system1_copy = deepcopy(rxn_system1)
    evolutionary_solve!(rxn_system1_copy, system1.K_eqs, CMAES(); reltol=1.0e-6)
    @test isapprox(rxn_system1_copy.concs, system1.concs_eq; rtol=1.0e-6)
end

@testset verbose=true "System 2" begin
    # Set up reaction system
    rxn_system2 = ReactionSystem(system2.stoich, copy(system2.concs), system2.K_eqs; φ=1.0)
    # Run and test NEKMC
    simulate!(
        rxn_system2;
        n_iter=Int(1e+9),
        chunk_iter=Int(1e+2),
        ε=1e-5,
        ε_scale=0.5,
        ε_concs=2.0,
        # tol_ε=0.0,
    )
    @test isapprox(rxn_system2.concs, system2.concs_eq; rtol=1.0e-4)
    # Run and test CMA solver
    rxn_system2_copy = deepcopy(rxn_system2)
    evolutionary_solve!(rxn_system2_copy, system2.K_eqs, CMAES(); maxiter=1000000, reltol=1.0e-6)
    @test isapprox(rxn_system2_copy.concs, system2.concs_eq; rtol=1.0e-6)
end
