using Printf
using VerifySDP

# Standard dimerized spinful Hubbard chain with weak inter-dimer hopping.
#
# This scan compares several selected moment families, including inter-dimer
# moments. It is used to document which moment choices have small response
# error and which show linear or quadratic small-theta mismatch behavior.

project_dir = normpath(joinpath(@__DIR__, ".."))
results_dir = joinpath(project_dir, "results")
mkpath(results_dir)

L = 4
U = 2.0
theta_grid = [-1e-2, -5e-3, 5e-3, 1e-2]
selected_sets = [
    :charge_imbalance,
    :bond_kinetic,
    :double_occupancy,
    :double_occupancy_sites,
    :double_occupancy_imbalance,
    :inter_bond_kinetic,
    :inter_density,
    :inter_charge_imbalance,
    :charge_and_bond,
]

function ratios_at(rows, abs_theta)
    return [row.ratio for row in rows if isapprox(abs(row.theta), abs_theta; atol = 1e-12)]
end

function min_or_nan(xs)
    isempty(xs) ? NaN : minimum(xs)
end

rows = NamedTuple[]

for selected in selected_sets
    @printf("START L=%d U=%.1f selected=%s\n", L, U, selected)
    try
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
        ratio_005 = min_or_nan(ratios_at(result.mismatch_rows, 5e-3))
        ratio_001 = min_or_nan(ratios_at(result.mismatch_rows, 1e-2))
        quadratic_005 = minimum(row.delta_norm / (row.theta * row.theta)
            for row in result.mismatch_rows if isapprox(abs(row.theta), 5e-3; atol = 1e-12))
        quadratic_001 = minimum(row.delta_norm / (row.theta * row.theta)
            for row in result.mismatch_rows if isapprox(abs(row.theta), 1e-2; atol = 1e-12))
        stability = ratio_005 / ratio_001
        statuses = join(unique(row.sdp_status for row in result.mismatch_rows), "|")
        row = (;
            L,
            U,
            selected = string(selected),
            moment_count = length(result.reference.sdp_moments),
            physical_gap = result.reference.physical_gap,
            exact_energy_signed = result.reference.energy_gap,
            exact_energy_abs = abs(result.reference.energy_gap),
            moment_error = result.reference.moment_error,
            response_diff_norm = result.response.diff_norm,
            ratio_abs_theta_005 = ratio_005,
            ratio_abs_theta_010 = ratio_001,
            ratio_stability_005_over_010 = stability,
            quadratic_abs_theta_005 = quadratic_005,
            quadratic_abs_theta_010 = quadratic_001,
            reference_status = result.reference.sdp_status,
            mismatch_statuses = statuses,
        )
        push!(rows, row)
        @printf("DONE  ratio005=%.6e ratio010=%.6e stability=%.3f response=%.3e moment=%.3e status=%s\n",
            row.ratio_abs_theta_005, row.ratio_abs_theta_010,
            row.ratio_stability_005_over_010, row.response_diff_norm,
            row.moment_error, row.reference_status)
    catch err
        row = (;
            L,
            U,
            selected = string(selected),
            moment_count = 0,
            physical_gap = NaN,
            exact_energy_signed = NaN,
            exact_energy_abs = NaN,
            moment_error = NaN,
            response_diff_norm = NaN,
            ratio_abs_theta_005 = NaN,
            ratio_abs_theta_010 = NaN,
            ratio_stability_005_over_010 = NaN,
            quadratic_abs_theta_005 = NaN,
            quadratic_abs_theta_010 = NaN,
            reference_status = "ERROR",
            mismatch_statuses = sprint(showerror, err),
        )
        push!(rows, row)
        @printf("ERROR selected=%s: %s\n", selected, sprint(showerror, err))
    end
end

csv_file = joinpath(results_dir, "standard_moments_mosek_scan.csv")
top_file = joinpath(results_dir, "standard_moments_mosek_scan_top.txt")

headers = [
    "L",
    "U",
    "selected",
    "moment_count",
    "physical_gap",
    "exact_energy_signed",
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

open(csv_file, "w") do io
    println(io, join(headers, ","))
    for row in rows
        println(io, join((getfield(row, Symbol(h)) for h in headers), ","))
    end
end

valid = filter(row ->
        isfinite(row.ratio_abs_theta_005) &&
        row.exact_energy_abs <= 1e-6 &&
        row.moment_error <= 1e-5 &&
        row.response_diff_norm <= 1e-5 &&
        row.reference_status == "OPTIMAL" &&
        row.mismatch_statuses == "OPTIMAL",
    rows)
sort!(valid, by = row -> row.ratio_abs_theta_005, rev = true)

open(top_file, "w") do io
    println(io, "Top MOSEK moment candidates for standard dimerized Hubbard")
    println(io, "L=4 U=2 theta=+/-0.005,+/-0.01")
    println(io, "stability near 1 suggests linear scaling; near 0.5 suggests quadratic scaling")
    for row in valid
        @printf(io,
            "selected=%s moments=%d ratio005=%.6e ratio010=%.6e stability=%.3f quad005=%.6e response=%.3e exact=%.3e moment=%.3e\n",
            row.selected, row.moment_count, row.ratio_abs_theta_005,
            row.ratio_abs_theta_010, row.ratio_stability_005_over_010,
            row.quadratic_abs_theta_005, row.response_diff_norm,
            row.exact_energy_abs, row.moment_error)
    end
end

println("Wrote ", csv_file)
println("Wrote ", top_file)
