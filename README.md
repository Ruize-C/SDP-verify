# SDP-verify

Julia code and numerical outputs for checking reference-model assumptions in a
reduced SDP relaxation for small fermionic lattice models.

The experiments here compare exact diagonalization against a reduced
two-site-marginal SDP. They were written to support the summary in
[`reports/teacher_summary_density_vs_standard.md`](reports/teacher_summary_density_vs_standard.md).

## What Is Checked

For a reference Hamiltonian \(H_0\), selected moments \(S\), and a weak
perturbation \(H(\theta)=H_0+\theta V\), the scripts check:

1. **Reference exactness:** whether the reduced SDP energy and selected moments
   agree with the physical ground state at \(H_0\).
2. **First-order response matching:** whether the SDP and physical moment maps
   have the same first derivative with respect to fields coupled to the
   selected moments.
3. **Mismatch lower bound:** whether
   \(\|\delta m_S(\theta)\|/|\theta|\) stays bounded away from zero for small
   \(\theta\).

## Models In This Repository

### Standard Dimerized Spinful Hubbard Chain

The reference \(H_0\) is a product of isolated two-site spinful Hubbard dimers.
The perturbation is weak inter-dimer hopping.

Main scripts:

- `scripts/scan_standard_moments_mosek.jl`
- `scripts/scan_standard_dimer_local_moments_mosek.jl`

### Density-Perturbed Extended Spinful Hubbard Dimer Chain

The reference \(H_0\) is the same isolated-dimer Hamiltonian. The perturbation
is an inter-dimer density-density coupling. The selected moments are the bond
kinetic observables on each strong dimer.

Main scripts:

- `scripts/check_density_hubbard_mosek_L4.jl`
- `scripts/check_density_hubbard_mosek_L6_light.jl`
- `scripts/check_density_hubbard_mosek_L6_response.jl`

## Repository Layout

- `src/VerifySDP.jl`: model definitions, exact diagonalization, moment
  extraction, and the reduced SDP builder.
- `scripts/`: reproducibility scripts for the two models discussed in the
  report.
- `reports/`: summary and detailed notes.
- `output/pdf/teacher_summary_density_vs_standard.pdf`: rendered PDF version
  of the teacher-facing summary.
- `results/`: numerical outputs used by the reports.
- `test/runtests.jl`: lightweight non-MOSEK smoke test.

## Setup

Install Julia, then instantiate the project:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

MOSEK is used for the reported high-accuracy SDP runs. To reproduce those
runs, install MOSEK/MosekTools and set `MOSEKLM_LICENSE_FILE` in your local
environment. Do not commit the license file.

For a quick license-free smoke test:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Reproducing The Reported Runs

From the repository root:

```bash
julia --project=. scripts/scan_standard_moments_mosek.jl
julia --project=. scripts/scan_standard_dimer_local_moments_mosek.jl
julia --project=. scripts/check_density_hubbard_mosek_L4.jl
julia --project=. scripts/check_density_hubbard_mosek_L6_light.jl
julia --project=. scripts/check_density_hubbard_mosek_L6_response.jl
```

The scripts write outputs into `results/`.

## Notes

- The systems are intentionally small because the goal is to compare exact
  ground-state data with SDP solutions.
- The code uses relative paths from each script location, so the repository can
  be moved without editing path constants.
- `Manifest.toml` is included for reproducibility of the Julia environment used
  during the numerical checks.
