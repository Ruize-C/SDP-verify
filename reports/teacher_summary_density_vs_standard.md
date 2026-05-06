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

The quantities in this subsection are diagnostics used in this numerical
report, not standard names from the Hubbard-model literature.

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

## Terminology And References

- The terms "Hubbard model", "spinful Hubbard chain", "hopping", "onsite
  interaction", "half filling", and "double occupancy" follow the standard
  Hubbard-model terminology; see Hubbard's original paper and standard
  treatments of the one-dimensional Hubbard model [[1](#references),
  [2](#references)].
- The words "dimerized" and "alternating" refer to chains with alternating
  strong and weak bonds; this terminology is standard in the dimerized or
  alternating Hubbard-chain literature [[3](#references), [4](#references)].
- The term "extended Hubbard model" refers here to adding an intersite
  density-density interaction to the onsite Hubbard interaction; this is the
  standard usage in the extended-Hubbard literature [[5](#references)].
- "Response matrix" is used in the usual static linear-response sense,
  following the Kubo linear-response terminology [[6](#references)].
- "Semidefinite programming" follows the standard convex-optimization
  terminology of Vandenberghe and Boyd [[7](#references)].

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

- Bond kinetic moments on strong dimers:
  \(M_r=K_{2r-1,2r}\).
- Double occupancy moments on strong dimers:
  \(M_r=D_{2r-1}+D_{2r}\), and also the site-resolved versions
  \(D_{2r-1}\), \(D_{2r}\).
- Spin correlation moments on strong dimers:
  \(S^z_{2r-1}S^z_{2r}\), and
  \(\mathbf S_{2r-1}\cdot\mathbf S_{2r}
  =S^z_{2r-1}S^z_{2r}
  +\frac12(S^+_{2r-1}S^-_{2r}+S^-_{2r-1}S^+_{2r})\).
- Onsite-pair hopping moments between the two sites of a strong dimer:
  \(P_r=c^\dagger_{2r-1,\uparrow}c^\dagger_{2r-1,\downarrow}
  c_{2r,\downarrow}c_{2r,\uparrow}+h.c.\).
- Singlet/triplet projectors on a strong dimer:
  \(\Pi_{s,r}=|s_r\rangle\langle s_r|\) and
  \(\Pi_{t0,r}=|t^0_r\rangle\langle t^0_r|\), where
  \(|s_r\rangle=(|\uparrow,\downarrow\rangle-|\downarrow,\uparrow\rangle)/\sqrt2\)
  and
  \(|t^0_r\rangle=(|\uparrow,\downarrow\rangle+|\downarrow,\uparrow\rangle)/\sqrt2\)
  on dimer \(D_r\).
- Bond kinetic moments on weak inter-dimer bonds:
  \(M_r=K_{2r,2r+1}\) on the weak bonds \(B_r\).

**Verification results.**

1. **Reference exactness.** The reference SDP was solved at \(H_0\). The
   following rows give representative moment choices from the scan.

   | Selected moments | Energy error | Moment error | Status |
   | --- | ---: | ---: | --- |
   | \(K_{2r-1,2r}\) | \(1.47\times10^{-9}\) | \(9.37\times10^{-8}\) | OPTIMAL |
   | \(D_{2r-1}+D_{2r}\) | \(1.47\times10^{-9}\) | \(4.74\times10^{-8}\) | OPTIMAL |
   | \(K_{2r,2r+1}\) | \(1.47\times10^{-9}\) | \(4.65\times10^{-16}\) | OPTIMAL |

2. **First-order response matching.** Strong-dimer moments had small
   response-matching errors. The weak-bond kinetic moment had a large
   response-matching error.

   | Selected moments | Support | Response error |
   | --- | --- | ---: |
   | \(K_{2r-1,2r}\) | strong dimers \(D_r\) | \(4.86\times10^{-7}\) |
   | \(D_{2r-1}+D_{2r}\) | strong dimers \(D_r\) | \(1.56\times10^{-7}\) |
   | \(\mathbf S_{2r-1}\cdot\mathbf S_{2r}\) | strong dimers \(D_r\) | \(1.44\times10^{-7}\) |
   | \(\Pi_{s,r}\) | strong dimers \(D_r\) | \(1.10\times10^{-7}\) |
   | \(K_{2r,2r+1}\) | weak bonds \(B_r\) | \(7.05\) |

3. **Mismatch lower bound.** For strong-dimer moments, the ratio
   \(\|\delta m_S(\theta)\|/|\theta|\) approximately doubled when
   \(|\theta|\) doubled from \(0.005\) to \(0.01\). This indicates quadratic,
   not linear, scaling in the tested range. The weak-bond kinetic moment had a
   stable linear ratio, but it failed the response-matching check above.

   | Selected moments | Support | Ratio at 0.005 | Ratio at 0.01 | Observation |
   | --- | --- | ---: | ---: | --- |
   | \(K_{2r-1,2r}\) | strong dimers | \(5.96\times10^{-3}\) | \(1.17\times10^{-2}\) | quadratic |
   | \(D_{2r-1}+D_{2r}\) | strong dimers | \(7.04\times10^{-3}\) | \(1.40\times10^{-2}\) | quadratic |
   | \(\mathbf S_{2r-1}\cdot\mathbf S_{2r}\) | strong dimers | \(6.36\times10^{-3}\) | \(1.26\times10^{-2}\) | quadratic |
   | \(P_r\) | strong dimers | \(5.42\times10^{-3}\) | \(1.09\times10^{-2}\) | quadratic |
   | \(K_{2r,2r+1}\) | weak bonds | \(4.59\) | \(4.59\) | linear ratio, response failed |

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

1. **Reference exactness.** At \(H_0\), both the SDP energy and selected
   bond-kinetic moments matched the exact ground-state values to high
   precision.

   | System | Energy error | Moment error | Status |
   | --- | ---: | ---: | --- |
   | \(L=4,U=4\) | \(5.74\times10^{-12}\) | \(1.01\times10^{-7}\) | OPTIMAL |
   | \(L=6,U=4\) | \(8.18\times10^{-13}\) | \(1.03\times10^{-9}\) | OPTIMAL |

2. **First-order response matching.** The SDP and physical response matrices
   were compared by finite differences with field step \(h=0.005\).

   | System | Selected moments | Response error |
   | --- | --- | ---: |
   | \(L=4,U=4\) | \(K_{2r-1,2r}\) | \(3.04\times10^{-5}\) |
   | \(L=6,U=4\) | \(K_{2r-1,2r}\) | \(3.46\times10^{-6}\) |

3. **Mismatch lower bound.** The ratio
   \(\|\delta m_S(\theta)\|/|\theta|\) stayed nearly constant over the tested
   \(\theta\)-range, giving a numerical lower-bound constant.

   | System | Tested \(|\theta|\) | Ratio range | Conservative \(c_{\mathrm{mis}}\) |
   | --- | --- | ---: | ---: |
   | \(L=4,U=4\) | \(0.002,0.005,0.01,0.02\) | \(0.124\)--\(0.125\) | \(0.12\) |
   | \(L=6,U=4\) | \(0.005,0.01\) | \(0.215\)--\(0.216\) | \(0.21\) |

## Observed Status

For the standard hopping perturbation, the reference exactness error was small,
the response check for moments supported on strong dimers passed numerically,
and the corresponding mismatch was quadratic in \(\theta\) in the tested cases.

For the density-density perturbation, the tested \(L=4\) and \(L=6\) cases both
gave small reference errors, small response-matching errors, and a linear
mismatch ratio bounded away from zero on the tested \(\theta\)-range.

## References

[1] J. Hubbard, "Electron correlations in narrow energy bands," Proceedings of
the Royal Society A 276, 238-257 (1963).
https://doi.org/10.1098/rspa.1963.0204

[2] F. H. L. Essler, H. Frahm, F. Göhmann, A. Klümper, and V. E. Korepin,
The One-Dimensional Hubbard Model, Cambridge University Press (2005).
https://doi.org/10.1017/CBO9780511534843

[3] K. Penc and F. Mila, "Charge gap in the one-dimensional dimerized Hubbard
model at quarter-filling," Physical Review B 50, 11429 (1994).
https://doi.org/10.1103/PhysRevB.50.11429

[4] S. R. White, R. M. Noack, and D. J. Scalapino, "Density-matrix
renormalization-group studies of the alternating Hubbard model," Physical
Review B 51, 10287 (1995). https://doi.org/10.1103/PhysRevB.51.10287

[5] E. Jeckelmann, "Ground-State Phase Diagram of a Half-Filled
One-Dimensional Extended Hubbard Model," Physical Review Letters 89, 236401
(2002). https://doi.org/10.1103/PhysRevLett.89.236401

[6] R. Kubo, "Statistical-Mechanical Theory of Irreversible Processes. I,"
Journal of the Physical Society of Japan 12, 570-586 (1957).
https://doi.org/10.1143/JPSJ.12.570

[7] L. Vandenberghe and S. Boyd, "Semidefinite Programming," SIAM Review 38,
49-95 (1996). https://doi.org/10.1137/1038003
