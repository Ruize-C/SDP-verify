using LinearAlgebra
using Printf
using VerifySDP

# Density-perturbed extended spinful Hubbard dimer chain, L=6.
#
# This script computes the full finite-difference response matrix at H0 and
# compares the SDP response with the exact physical response.

function print_vector(io, label, v)
    println(io, label, " = [", join((@sprintf("%.10e", x) for x in v), ", "), "]")
end

function print_matrix(io, label, A)
    println(io, label, " =")
    for i in axes(A, 1)
        println(io, "  [", join((@sprintf("%.10e", A[i, j]) for j in axes(A, 2)), ", "), "]")
    end
end

function finite_response(f, s::Int; eps::Float64)
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
    L = 6,
    U = 4.0,
    perturbation = :density,
    selected = :bond_kinetic,
)

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

summary_file = joinpath(results_dir, "density_hubbard_L6_U4_bond_mosek_response_summary.txt")

open(summary_file, "w") do io
    println(io, "model=weakly_coupled_extended_spinful_hubbard_dimer_chain")
    println(io, "L=6")
    println(io, "U=4.0")
    println(io, "perturbation=density")
    println(io, "selected=bond_kinetic")
    println(io, "solver=MOSEK")
    println(io, "reference_solver_eps_abs=1e-10")
    println(io, "reference_solver_eps_rel=1e-10")
    println(io, "response_field_step=", response_eps_h)
    println(io, "response_solver_eps_abs=1e-10")
    println(io, "response_solver_eps_rel=1e-10")
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

println("== Density Hubbard L=6 MOSEK response check ==")
@printf("reference sdp - physical  %.12e\n", sdp0.energy - phys0.energy)
@printf("moment error norm          %.12e\n", norm(sdp0.moments .- phys0.moments))
@printf("response diff norm         %.12e\n", norm(sdp_response - phys_response))
println("Wrote ", summary_file)
