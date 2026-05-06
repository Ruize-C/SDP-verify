using LinearAlgebra
using Printf
using VerifySDP

# Density-perturbed extended spinful Hubbard dimer chain, L=6.
#
# This script checks reference exactness and small-theta mismatch ratios.
# The full L=6 response matrix is computed by
# `check_density_hubbard_mosek_L6_response.jl`.

project_dir = normpath(joinpath(@__DIR__, ".."))
results_dir = joinpath(project_dir, "results")
mkpath(results_dir)

model = spinful_hubbard_dimer_model(
    L = 6,
    U = 4.0,
    perturbation = :density,
    selected = :bond_kinetic,
)

theta_grid = [-1e-2, -5e-3, 5e-3, 1e-2]
summary_file = joinpath(results_dir, "density_hubbard_L6_U4_bond_mosek_mismatch_summary.txt")
mismatch_file = joinpath(results_dir, "density_hubbard_L6_U4_bond_mosek_mismatch.csv")

phys0 = physical_observables(model)
sdp0 = solve_sdp_observables(model; solver = :mosek, eps_abs = 1e-8, eps_rel = 1e-8, max_iters = 300)
ref_moment_error = norm(sdp0.moments .- phys0.moments)

open(summary_file, "w") do io
    println(io, "model=weakly_coupled_extended_spinful_hubbard_dimer_chain")
    println(io, "L=6")
    println(io, "U=4.0")
    println(io, "perturbation=density")
    println(io, "selected=bond_kinetic")
    println(io, "solver=MOSEK eps_abs=eps_rel=1e-8 max_iters=300")
    println(io, "response=not_run")
    println(io, "reference_phys_energy=", phys0.energy)
    println(io, "reference_sdp_energy=", sdp0.energy)
    println(io, "reference_sdp_minus_phys=", sdp0.energy - phys0.energy)
    println(io, "reference_physical_gap=", phys0.gap)
    println(io, "reference_moment_error_norm=", ref_moment_error)
    println(io, "reference_sdp_status=", sdp0.termination)
end

open(mismatch_file, "w") do io
    println(io, "theta,abs_theta,delta_norm,delta_norm_over_abs_theta,delta_norm_over_theta_squared,sdp_energy,physical_energy,sdp_minus_physical,sdp_status")
    for theta in theta_grid
        phys = physical_observables(model; theta)
        sdp = solve_sdp_observables(model; theta, solver = :mosek, eps_abs = 1e-8, eps_rel = 1e-8, max_iters = 300)
        delta = sdp.moments .- phys.moments
        delta_norm = norm(delta)
        println(io, join((
            theta,
            abs(theta),
            delta_norm,
            delta_norm / abs(theta),
            delta_norm / (theta * theta),
            sdp.energy,
            phys.energy,
            sdp.energy - phys.energy,
            sdp.termination,
        ), ","))
    end
end

println("== Density Hubbard L=6 MOSEK mismatch check ==")
@printf("reference sdp - physical  %.12e\n", sdp0.energy - phys0.energy)
@printf("moment error norm          %.12e\n", ref_moment_error)
println("theta,delta_norm,ratio,quadratic_ratio,status")
for line in eachline(mismatch_file)
    startswith(line, "theta,") && continue
    parts = split(line, ",")
    @printf("%s,%s,%s,%s,%s\n", parts[1], parts[3], parts[4], parts[5], parts[9])
end
println("Wrote ", summary_file)
println("Wrote ", mismatch_file)
