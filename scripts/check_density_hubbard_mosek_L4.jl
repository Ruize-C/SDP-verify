using LinearAlgebra
using Printf
using VerifySDP

# Density-perturbed extended spinful Hubbard dimer chain, L=4.
#
# This script produces the L=4 MOSEK numbers used in the summary:
#   * reference exactness at H0;
#   * first-order response matching for dimer bond-kinetic moments;
#   * small-theta mismatch ratios for the density-density perturbation.

function print_vector(io, label, v)
    println(io, label, " = [", join((@sprintf("%.10e", x) for x in v), ", "), "]")
end

function print_matrix(io, label, A)
    println(io, label, " =")
    for i in axes(A, 1)
        println(io, "  [", join((@sprintf("%.10e", A[i, j]) for j in axes(A, 2)), ", "), "]")
    end
end

function finite_response(f, s::Int; eps::Float64 = 1e-4)
    out0 = f(zeros(s))
    chi = zeros(length(out0), s)
    for a in 1:s
        hp = zeros(s)
        hm = zeros(s)
        hp[a] = eps
        hm[a] = -eps
        chi[:, a] = (f(hp) - f(hm)) ./ (2eps)
    end
    return chi
end

project_dir = normpath(joinpath(@__DIR__, ".."))
results_dir = joinpath(project_dir, "results")
mkpath(results_dir)

model = spinful_hubbard_dimer_model(
    L = 4,
    U = 4.0,
    perturbation = :density,
    selected = :bond_kinetic,
)

theta_grid = [-2e-2, -1e-2, -5e-3, -2e-3, 2e-3, 5e-3, 1e-2, 2e-2]
response_eps_h = 5e-3
s = length(model.selected)

phys0 = physical_observables(model)
sdp0 = solve_sdp_observables(model; solver = :mosek, eps_abs = 1e-10, eps_rel = 1e-10, max_iters = 1000)

phys_response = finite_response(
    h -> physical_observables(model; theta = 0.0, fields = h).moments,
    s,
    eps = response_eps_h,
)
sdp_response = finite_response(
    h -> solve_sdp_observables(model; theta = 0.0, fields = h,
        solver = :mosek, eps_abs = 1e-10, eps_rel = 1e-10, max_iters = 1000).moments,
    s,
    eps = response_eps_h,
)

summary_file = joinpath(results_dir, "density_hubbard_L4_U4_bond_mosek_summary.txt")
mismatch_file = joinpath(results_dir, "density_hubbard_L4_U4_bond_mosek_mismatch.csv")
components_file = joinpath(results_dir, "density_hubbard_L4_U4_bond_mosek_delta_components.csv")
coeff_file = joinpath(results_dir, "density_hubbard_L4_U4_bond_mosek_leading_coefficients.csv")

theta_rows = NamedTuple[]
for theta in theta_grid
    phys = physical_observables(model; theta)
    sdp = solve_sdp_observables(model; theta, solver = :mosek, eps_abs = 1e-8, eps_rel = 1e-8, max_iters = 200)
    delta = sdp.moments .- phys.moments
    push!(theta_rows, (;
        theta,
        phys_energy = phys.energy,
        sdp_energy = sdp.energy,
        phys_moments = phys.moments,
        sdp_moments = sdp.moments,
        delta,
        delta_norm = norm(delta),
        ratio = norm(delta) / abs(theta),
        quadratic_ratio = norm(delta) / (theta * theta),
        status = sdp.termination,
    ))
end

open(summary_file, "w") do io
    println(io, "model=weakly_coupled_extended_spinful_hubbard_dimer_chain")
    println(io, "L=4")
    println(io, "U=4.0")
    println(io, "perturbation=density")
    println(io, "selected=bond_kinetic")
    println(io, "solver=MOSEK")
    println(io, "reference_solver_eps_abs=1e-10")
    println(io, "reference_solver_eps_rel=1e-10")
    println(io, "response_field_step=", response_eps_h)
    println(io, "response_solver_eps_abs=1e-10")
    println(io, "response_solver_eps_rel=1e-10")
    println(io, "mismatch_solver_eps_abs=1e-8")
    println(io, "mismatch_solver_eps_rel=1e-8")
    println(io, "reference_phys_energy=", phys0.energy)
    println(io, "reference_sdp_energy=", sdp0.energy)
    println(io, "reference_sdp_minus_phys=", sdp0.energy - phys0.energy)
    println(io, "reference_physical_gap=", phys0.gap)
    println(io, "reference_moment_error_norm=", norm(sdp0.moments .- phys0.moments))
    println(io, "reference_sdp_status=", sdp0.termination)
    println(io, "response_diff_norm=", norm(sdp_response - phys_response))
    print_vector(io, "physical_moments", phys0.moments)
    print_vector(io, "sdp_moments", sdp0.moments)
    print_matrix(io, "physical_chi", phys_response)
    print_matrix(io, "sdp_chi", sdp_response)
end

open(mismatch_file, "w") do io
    println(io, "theta,abs_theta,delta_norm,delta_norm_over_abs_theta,delta_norm_over_theta_squared,sdp_energy,physical_energy,sdp_minus_physical,sdp_status")
    for row in theta_rows
        println(io, join((
            row.theta,
            abs(row.theta),
            row.delta_norm,
            row.ratio,
            row.quadratic_ratio,
            row.sdp_energy,
            row.phys_energy,
            row.sdp_energy - row.phys_energy,
            row.status,
        ), ","))
    end
end

open(components_file, "w") do io
    headers = ["theta"]
    for spec in model.selected
        push!(headers, "phys_" * spec.name)
        push!(headers, "sdp_" * spec.name)
        push!(headers, "delta_" * spec.name)
        push!(headers, "delta_over_theta_" * spec.name)
    end
    println(io, join(headers, ","))
    for row in theta_rows
        values = Any[row.theta]
        for k in eachindex(model.selected)
            delta = row.sdp_moments[k] - row.phys_moments[k]
            append!(values, (row.phys_moments[k], row.sdp_moments[k], delta, delta / row.theta))
        end
        println(io, join(values, ","))
    end
end

open(coeff_file, "w") do io
    println(io, "abs_theta,central_delta_derivative_norm,central_delta_derivative_components")
    for h in [2e-3, 5e-3, 1e-2, 2e-2]
        pos = only(row for row in theta_rows if isapprox(row.theta, h; atol = 1e-14))
        neg = only(row for row in theta_rows if isapprox(row.theta, -h; atol = 1e-14))
        coeff = (pos.delta .- neg.delta) ./ (2h)
        println(io, join((h, norm(coeff), "[" * join(coeff, ";") * "]"), ","))
    end
end

println("== Density Hubbard L=4 MOSEK check ==")
@printf("reference sdp - physical  %.12e\n", sdp0.energy - phys0.energy)
@printf("moment error norm          %.12e\n", norm(sdp0.moments .- phys0.moments))
@printf("response diff norm         %.12e\n", norm(sdp_response - phys_response))
println("theta,delta_norm,ratio,status")
for row in theta_rows
    @printf("%.8e,%.12e,%.12e,%s\n", row.theta, row.delta_norm, row.ratio, row.status)
end
println("Wrote ", summary_file)
println("Wrote ", mismatch_file)
println("Wrote ", components_file)
println("Wrote ", coeff_file)
