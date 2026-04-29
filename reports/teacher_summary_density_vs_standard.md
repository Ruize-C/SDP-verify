# Summary of Current Model Verification Study

## Goal

The goal of this round was to find a concrete fermionic model that can support
the perturbative SDP improvement framework in the draft. In particular, we
wanted to verify three points:

1. exactness of the reduced SDP at a reference Hamiltonian \(H_0\);
2. first-order response matching between the reduced SDP and the physical
   ground state at \(H_0\);
3. a nontrivial lower bound
   \[
   \|\delta m_S(\theta)\| \ge c_{\mathrm{mis}} |\theta|
   \]
   for selected moments after a weak perturbation
   \(H(\theta)=H_0+\theta V\).

We tested two spinful Hubbard-type candidates, both based on an isolated
Hubbard-dimer reference point, i.e. a decoupled product of two-site spinful
Hubbard dimers. Here \(L\) denotes the number of lattice sites and \(U\) is the
onsite Hubbard interaction.

## Reported Quantities

- Reference energy error:
  \(|E_{\mathrm{SDP}}(H_0)-E_{\mathrm{phys}}(H_0)|\).
- Selected-moment error:
  \(\|m_S^{\mathrm{SDP}}(0)-m_S^{\mathrm{phys}}(0)\|\).
- Response-matching error:
  \(\|D_hm_S^{\mathrm{SDP}}(0)-D_hm_S^{\mathrm{phys}}(0)\|\), the norm of the
  difference between the SDP and physical first-order response matrices.
- Mismatch ratio:
  \(\|\delta m_S(\theta)\|/|\theta|\), where
  \(\delta m_S(\theta)=m_S^{\mathrm{SDP}}(\theta)-m_S^{\mathrm{phys}}(\theta)\).
  This ratio is used to estimate \(c_{\mathrm{mis}}\).
- Solver:
  all reported SDP values below were computed with MOSEK.

## Model 1: Standard Dimerized Spinful Hubbard Chain

**Model.** The first model is the standard alternating-hopping spinful Hubbard
chain. The reference Hamiltonian \(H_0\) is a product of isolated Hubbard
dimers, and the weak perturbation is inter-dimer hopping:

\[
V_{\mathrm{hop}}
=-\sum_{k,\sigma}
(c^\dagger_{2k,\sigma}c_{2k+1,\sigma}+h.c.).
\]

Here "dimer-local moments" refers to observables supported inside each strong
dimer, i.e. inside each isolated two-site unit of \(H_0\), such as bond kinetic
energy, double occupancy, and spin correlations.
"Inter-dimer bond moments" refers to hopping-type observables on the weak
bonds connecting neighboring dimers.

**Verification results.**

1. **Reference exactness.** With MOSEK, the reference energy error was around
   \(10^{-9}\).

2. **First-order response matching.** Dimer-local selected moments matched well
   at \(H_0\). Inter-dimer bond moments gave a large response-matching error.

3. **Mismatch lower bound.** For all dimer-local observables tested, including
   bond kinetic, double occupancy, spin correlations, pair hopping, and
   singlet/triplet projectors, the mismatch scaled quadratically:

\[
\|\delta m_S(\theta)\| \sim C\theta^2,
\]

   rather than linearly in \(|\theta|\). Inter-dimer bond moments did show a
   strong linear mismatch, while response matching failed for those moments.

## Model 2: Density-Perturbed Extended Spinful Hubbard Dimer Chain

**Model.** The second model uses the same isolated Hubbard-dimer reference
\(H_0\), but replaces weak hopping by an inter-dimer density-density
perturbation:

\[
V_{\mathrm{dens}}
=\sum_k n_{2k}n_{2k+1}.
\]

Here \(n_i=n_{i,\uparrow}+n_{i,\downarrow}\) is the total occupation at site
\(i\).

The selected moments are the bond kinetic observables on each strong dimer:

\[
M_k=\sum_\sigma
(c^\dagger_{2k-1,\sigma}c_{2k,\sigma}+h.c.).
\]

**Verification results.**

1. **Reference exactness.** With MOSEK, the reference point was exact to high
   precision.
   - \(L=4, U=4\): reference energy error about \(5.7\times10^{-12}\);
     selected-moment error about \(1.0\times10^{-7}\).
   - \(L=6, U=4\): reference energy error about \(8.2\times10^{-13}\);
     selected-moment error about \(1.0\times10^{-9}\).

2. **First-order response matching.** The SDP and physical response matrices
   had small response-matching errors.
   - \(L=4, U=4\): response-matching error about \(3.0\times10^{-5}\).
   - \(L=6, U=4\): response-matching error about \(3.5\times10^{-6}\).

3. **Mismatch lower bound.** The selected bond-kinetic moments showed a stable
   linear mismatch after the density perturbation.
   - \(L=4, U=4\): \(\|\delta m_S(\theta)\|/|\theta|\approx 0.125\) for
     \(|\theta|\le 0.02\), corresponding to
     \(c_{\mathrm{mis}}\approx 0.12\) on the tested range.
   - \(L=6, U=4\): \(\|\delta m_S(\theta)\|/|\theta|\approx
     0.215\text{--}0.216\) for \(|\theta|=0.005,0.01\), corresponding to
     \(c_{\mathrm{mis}}\approx 0.21\) on the tested range.

## Observed Status

For the standard hopping perturbation, the reference exactness error was small,
the dimer-local response check passed numerically, and the dimer-local mismatch
was quadratic in \(\theta\) in the tested cases.

For the density-density perturbation, the tested \(L=4\) and \(L=6\) cases both
gave small reference errors, small response-matching errors, and a linear
mismatch ratio bounded away from zero on the tested \(\theta\)-range.
