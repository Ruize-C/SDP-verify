using Printf
using VerifySDP

# Standard dimerized spinful Hubbard chain with weak inter-dimer hopping.
#
# This scan restricts attention to dimer-local selected moments, i.e. moments
# supported inside the strong dimers of the reference Hamiltonian H0.

project_dir = normpath(joinpath(@__DIR__, ".."))
results_dir = joinpath(project_dir, "results")
mkpath(results_dir)

L = 4
U = 2.0
theta_grid = [-1e-2, -5e-3, 5e-3, 1e-2]

selected_sets = [
    :bond_kinetic,
    :double_occupancy,
    :double_occupancy_sites,
    :double_occupancy_imbalance,
    :dimer_spin_z_corr,
    :dimer_spin_exchange,
    :dimer_spin_dot,
    :dimer_pair_hopping,
    :dimer_singlet_projector,
    :dimer_triplet0_projector,
]

function min_ratio(rows, abs_theta)
    values = [row.ratio for row in rows if isapprox(abs(row.theta), abs_theta; atol = 1e-12)]
    return isempty(values) ? NaN : minimum(values)
end

function min_quad(rows, abs_theta)
    values = [row.delta_norm / (row.theta * row.theta)
        for row in rows if isapprox(abs(row.theta), abs_theta; atol = 1e-12)]
    return isempty(values) ? NaN : minimum(values)
end

rows = NamedTuple[]

for selected in selected_sets
    @printf("START selected=%s\n", selected)
    result = run_spinful_hubbard_probe(;
        L,
        U,
        perturbation = :hopping,
        selected,
        theta_grid = theta_grid,
        solver = :mosek,
        eps_abs = 1e-8,
        eps_rel = 1e-8,
        max_iters = 200,
    )
    ratio005 = min_ratio(result.mismatch_rows, 5e-3)
    ratio010 = min_ratio(result.mismatch_rows, 1e-2)
    row = (;
        L,
        U,
        selected = string(selected),
        moment_count = length(result.reference.sdp_moments),
        exact_energy_abs = abs(result.reference.energy_gap),
        moment_error = result.reference.moment_error,
        response_diff_norm = result.response.diff_norm,
        ratio_abs_theta_005 = ratio005,
        ratio_abs_theta_010 = ratio010,
        ratio_stability_005_over_010 = ratio005 / ratio010,
        quadratic_abs_theta_005 = min_quad(result.mismatch_rows, 5e-3),
        quadratic_abs_theta_010 = min_quad(result.mismatch_rows, 1e-2),
        reference_status = result.reference.sdp_status,
        mismatch_statuses = join(unique(row.sdp_status for row in result.mismatch_rows), "|"),
    )
    push!(rows, row)
    @printf("DONE selected=%s ratio005=%.6e ratio010=%.6e stability=%.3f response=%.3e moment=%.3e\n",
        row.selected, row.ratio_abs_theta_005, row.ratio_abs_theta_010,
        row.ratio_stability_005_over_010, row.response_diff_norm, row.moment_error)
end

headers = [
    "L",
    "U",
    "selected",
    "moment_count",
    "exact_energy_abs",
    "moment_error",
    "response_diff_norm",
    "ratio_abs_theta_005",
    "ratio_abs_theta_010",
    "ratio_stability_005_over_010",
    "quadratic_abs_theta_005",
    "quadratic_abs_theta_010",
    "reference_status",
    "mismatch_statuses",
]

csv_file = joinpath(results_dir, "standard_dimer_local_moments_mosek_scan.csv")
top_file = joinpath(results_dir, "standard_dimer_local_moments_mosek_scan_top.txt")

open(csv_file, "w") do io
    println(io, join(headers, ","))
    for row in rows
        println(io, join((getfield(row, Symbol(h)) for h in headers), ","))
    end
end

valid = filter(row ->
        row.exact_energy_abs <= 1e-6 &&
        row.moment_error <= 1e-5 &&
        row.response_diff_norm <= 1e-5 &&
        row.reference_status == "OPTIMAL" &&
        row.mismatch_statuses == "OPTIMAL",
    rows)
sort!(valid, by = row -> row.ratio_stability_005_over_010, rev = true)

open(top_file, "w") do io
    println(io, "Dimer-local MOSEK moment scan for standard dimerized Hubbard")
    println(io, "L=4 U=2 theta=+/-0.005,+/-0.01")
    println(io, "stability near 1 suggests first-order; near 0.5 suggests second-order")
    for row in valid
        @printf(io,
            "selected=%s ratio005=%.6e ratio010=%.6e stability=%.3f quad005=%.6e response=%.3e exact=%.3e moment=%.3e\n",
            row.selected, row.ratio_abs_theta_005, row.ratio_abs_theta_010,
            row.ratio_stability_005_over_010, row.quadratic_abs_theta_005,
            row.response_diff_norm, row.exact_energy_abs, row.moment_error)
    end
end

println("Wrote ", csv_file)
println("Wrote ", top_file)
