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
Hubbard-dimer reference point.

## Common Setup And Notation

- Lattice: an open one-dimensional chain with \(L\) sites, \(L\) even.
- Filling: half filling, \(N_e=L\) electrons.
- Strong dimers: \(D_r=(2r-1,2r)\), \(r=1,\dots,L/2\).
- Weak inter-dimer bonds: \(B_r=(2r,2r+1)\), \(r=1,\dots,L/2-1\).
- Operators:
  - \(c_{i,\sigma},c^\dagger_{i,\sigma}\): fermion annihilation/creation
    operators at site \(i\) with spin \(\sigma\in\{\uparrow,\downarrow\}\).
  - \(n_{i,\sigma}=c^\dagger_{i,\sigma}c_{i,\sigma}\), and
    \(n_i=n_{i,\uparrow}+n_{i,\downarrow}\).
  - \(D_i=n_{i,\uparrow}n_{i,\downarrow}\), the onsite double occupancy.
  - \(K_{ij}=\sum_\sigma(c^\dagger_{i,\sigma}c_{j,\sigma}
    +c^\dagger_{j,\sigma}c_{i,\sigma})\), the bond kinetic observable.
  - \(S_i^z=(n_{i,\uparrow}-n_{i,\downarrow})/2\),
    \(S_i^+=c^\dagger_{i,\uparrow}c_{i,\downarrow}\), and
    \(S_i^-=c^\dagger_{i,\downarrow}c_{i,\uparrow}\).
  - \(h.c.\) denotes the Hermitian conjugate.

The common reference Hamiltonian is the isolated-dimer spinful Hubbard model
with \(t_d=1\):

\[
H_0=\sum_{r=1}^{L/2}
\left[-K_{2r-1,2r}
+U(D_{2r-1}+D_{2r})\right].
\]

For response matching, fields are coupled to the selected moments
\(M_1,\dots,M_s\) by
\[
H(\theta,h)=H_0+\theta V-\sum_{a=1}^s h_aM_a.
\]

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
chain. It uses the common reference Hamiltonian \(H_0\) above and the weak
inter-dimer hopping perturbation

\[
V_{\mathrm{hop}}
=-\sum_{r=1}^{L/2-1}K_{2r,2r+1}.
\]

The reported scan used \(L=4\), \(U=2\), \(t_d=1\), and
\(\theta=\pm0.005,\pm0.01\).

The selected moment families tested were:

- Dimer-local bond kinetic moments:
  \(M_r=K_{2r-1,2r}\).
- Dimer-local double occupancy moments:
  \(M_r=D_{2r-1}+D_{2r}\), and also the site-resolved versions
  \(D_{2r-1}\), \(D_{2r}\).
- Dimer-local spin moments:
  \(S^z_{2r-1}S^z_{2r}\), and
  \(\mathbf S_{2r-1}\cdot\mathbf S_{2r}
  =S^z_{2r-1}S^z_{2r}
  +\frac12(S^+_{2r-1}S^-_{2r}+S^-_{2r-1}S^+_{2r})\).
- Dimer-local pair-hopping moments:
  \(P_r=c^\dagger_{2r-1,\uparrow}c^\dagger_{2r-1,\downarrow}
  c_{2r,\downarrow}c_{2r,\uparrow}+h.c.\).
- Dimer-local singlet/triplet projectors:
  \(\Pi_{s,r}=|s_r\rangle\langle s_r|\) and
  \(\Pi_{t0,r}=|t^0_r\rangle\langle t^0_r|\), where
  \(|s_r\rangle=(|\uparrow,\downarrow\rangle-|\downarrow,\uparrow\rangle)/\sqrt2\)
  and
  \(|t^0_r\rangle=(|\uparrow,\downarrow\rangle+|\downarrow,\uparrow\rangle)/\sqrt2\)
  on dimer \(D_r\).
- Inter-dimer bond kinetic moments:
  \(M_r=K_{2r,2r+1}\) on the weak bonds \(B_r\).

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
=\sum_{r=1}^{L/2-1} n_{2r}n_{2r+1}.
\]

The selected moments are the strong-dimer bond kinetic observables:

\[
M_r=K_{2r-1,2r},\qquad r=1,\dots,L/2.
\]

The reported checks used \(U=4\), \(t_d=1\), and \(L=4,6\).

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
