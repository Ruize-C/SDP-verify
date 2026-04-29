using Test
using LinearAlgebra
using VerifySDP

@testset "spinless dimer smoke test" begin
    model = dimer_model(L = 4, perturbation = :density, selected = :dimer_hopping)

    phys = physical_observables(model)
    @test isfinite(phys.energy)
    @test length(phys.moments) == 2
    @test phys.gap > 0

    sdp = solve_sdp_observables(model; solver = :scs, eps_abs = 1e-5, eps_rel = 1e-5, max_iters = 5_000)
    @test isfinite(sdp.energy)
    @test length(sdp.moments) == length(phys.moments)
end

@testset "spinful Hubbard model construction" begin
    model = spinful_hubbard_dimer_model(
        L = 4,
        U = 4.0,
        perturbation = :density,
        selected = :bond_kinetic,
    )

    phys = physical_observables(model)
    @test isfinite(phys.energy)
    @test length(phys.moments) == 2
    @test phys.gap > 0
end
