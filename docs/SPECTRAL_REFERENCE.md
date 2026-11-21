# Analytic spectral-reference interface

Every public matrix generator in the unitary Anymatrix group accepts an optional
terminal 'spectral' mode. It returns a MATLAB structure with known reference
spectral information **without constructing the matrix**.

## Example

~~~matlab
U = anymatrix('unitary/dirac_automaton',64,.45);
op = anymatrix('unitary/dirac_automaton',64,.45,'operator');
S = anymatrix('unitary/dirac_automaton',64,.45,'spectral');

theta = linspace(-pi,pi,400);
rho = S.limit.density(theta);         % analytic limiting trace density
lambda = S.finite.eigenvalues();      % exact 128 phases, computed on demand
moments = S.finite.moments(0:6);      % normalized finite trace moments
~~~

Use the **complete positional generator signature**, supplying empty
optional arguments as [] where necessary. The new mode occupies the same
terminal argument slot used by 'operator' in generators that already have it.
For generators without 'operator', it is an additional terminal argument.

~~~matlab
S = anymatrix('unitary/jacobi_szego',32,0,0,'cmv',1,[],'spectral');
S = anymatrix('unitary/orthogonal_blocks',21,[],[],'spectral');
S = anymatrix('unitary/haar_unitary',100,[],'spectral');
~~~

Direct MATLAB generator calls also work. The existing matrix output, two-output
properties interface, operator mode, and Anymatrix property scanning are unchanged.

## Structure contract

| Field | Meaning |
|---|---|
| S.family | Generator name. |
| S.dimension | Matrix order, not necessarily physical site count. |
| S.finite | Equal-weight finite-dimensional **trace spectral measure**. |
| S.limit | Large-order limiting **trace spectral measure** if proved or available. |
| S.ensemble | Expected spectral law of a random-matrix **ensemble**, not a realized matrix. |
| S.underlying | **Infinite-dimensional OPUC distinguished-vector measure**, where available; not the finite CMV trace law. |
| S.reference | Mathematical provenance where supplied. |
| S.notes | Qualification of known and unknown information. |

Finite data provide:
- finite.kind ('exact' or 'unknown');
- finite.eigenvalues: zero-argument function handle returning analytically known
  finite eigenvalues when available; an empty array otherwise;
- finite.moments: handle accepting integer powers k, returning
  tr(U^k)/dimension (shape follows k);
- finite.atoms: structure with fields angles and normalized weights, where
  multiplicity data are explicitly recorded.

The limit, ensemble, and underlying fields use a common structure:
kind, density, support, atoms, moments. Unknown values are [].
Densities, when provided, are with respect to dtheta on [-pi,pi], so
their integral over a period is 1. A 'support_only' descriptor does
**not** imply an analytic phase density. An atomic law has atoms.angles
and atoms.weights instead of a Lebesgue density.

## Implemented reference information

| Family | Finite reference | Other reference |
|---|---|---|
| cycle_shift and weighted_shift | Exact cyclic phases | Uniform large-size counting measure |
| diagonal_phase and phase_circulant | Prescribed finite phases | Uniform only for default equally spaced phases |
| orthogonal_blocks | Exact rotation pairs and possible +/-1 | No universal limit |
| cycle_propagator | Exact Fourier-mode phases | Wrapped arcsine density; Bessel moments |
| dirac_automaton | Exact two-branch phases | Two-band density (0<m<1), uniform m=0, atomic m=1 |
| grover_operator | Exact atoms and multiplicities | No universal continuum density |
| discrete_fourier | Exact n-dependent +/-1,+/-i multiplicities | Atomic, not a uniform phase law |
| walsh_hadamard, natural order | Exact +/-1 spectrum | Sequency permutation not generally involutory |
| time_frequency_shift | Exact monomial-cycle eigenvalues | No universal limit |
| jaynes_cummings | Exact excitation doublets plus cutoff boundary levels | No universal limit |
| haar_unitary and circular_beta | Sample spectrum not prescribed | Uniform expected density for finite ensembles |
| haar_orthogonal | Sample spectrum not prescribed | Asymptotically uniform counting measure |
| jacobi_szego, bernstein_szego, circular_jacobi, rogers_szego | No general finite counting formula | Known **underlying** circle weight |
| geronimus_cmv and period-one periodic_cmv | Boundary-dependent finite spectrum | Essential support bands; not an explicit density |

All other generators return a classified descriptor with unknown information
rather than fabricating spectral formulas or silently calling eig().

## Mathematical distinctions

For a unitary $U\in\mathbb C^{N\times N}$, the finite empirical measure is
written using uppercase $N$, the actual matrix order (not necessarily a
physical parameter named $n$):

```math
 \mu_N=\frac{1}{N}\sum_{j=1}^{N}\delta_{\arg\lambda_j},
 \qquad m_k=\frac{\mathop{\mathrm{tr}}(U^k)}{N}.
```

It is atomic regardless of the presence of a smooth limiting density.
For quantum walks and OPUC models, limit and finite quantities may differ,
and finite boundary phases can place eigenvalues in gaps of an infinite
periodic operator. For random matrices the ensemble mean is a different
object from a sample's histogram.

OPUC reference weights such as the Jacobi--Szegő lift describe the infinite
multiplication operator at its distinguished vector. A finite CMV
realization has a discrete Gauss--Szegő vector quadrature measure, generally
with unequal weights, while its empirical trace measure assigns weight 1/N
to each eigenvalue. These three objects must not be conflated.

The Dirac automaton and wrapped arcsine densities have square-root
singularities at band edges. Numerical integration should split intervals
at those endpoints rather than evaluate a trapezoidal rule on a grid
containing infinite endpoint values.

## Regression testing

The MATLAB suite explicitly enumerates 87 spectral-reference cases (other registered generators have separate tests),
analytically known finite eigenvalues against independent eig() results,
Fourier moments, DFT and Grover multiplicities, Dirac band density and
endpoint cases, wrapped phase periodicity, separation of OPUC measures,
random-number-state preservation, deferred large-order evaluation,
and unmodified legacy output modes.

The GitHub Actions workflow also executes the suite in an installed
Anymatrix external group. If CI has not completed successfully, do not
interpret the presence of tests as evidence that they passed.

**Source:** unitary_spectral_reference.m; this is an extension within the
unitary group, not a patch to upstream Higham--Mikaitis Anymatrix.

For `geronimus_cmv`, `S.underlying.density(theta)` is the **infinite**
constant-Verblunsky distinguished-vector absolutely continuous density.
`S.underlying.atoms` records an additional atom when present, and
`S.underlying.poisson_density(theta,r)` analytically smooths their sum
for `0<=r<1`. None of these is the finite equal-weight trace measure;
finite CMV eigenphases and cyclic weights must be computed separately.

`spherical_slepian` reports an unknown orthogonal-matrix trace phase law;
concentration eigenvalues belong instead to the distinct symmetric cap
matrix returned in `concentration` mode.

`fourier_laguerre_ball(P,L,tau,'spectral')` reports unknown finite trace
eigenphases for the P*L square radial-zonal orthogonal matrix. The Gaussian
radial and polar quadrature nodes do NOT represent matrix eigenvalues;
`'sampling'` gives a separate rectangular full-spherical isometry.

`abl_ladik_cmv` exposes fixed *finite trace* eigenphases and moments via
`S.finite` and varying finite *cyclic-vector* spectral weights in `S.cyclic`.
For each harmonic k, `S.flow.generator_values` describes the exact
Nenciu/Gekhtman-Nenciu isospectral weight evolution; no continuous
finite spectral density is inferred from the atomic measure.

`etna_hessenberg(n,'7.1',...)` provides no closed-form trace eigenphase
law; its common core parameters are sampled once per realization.
`etna_hessenberg(n,'7.5',...)` is a **modern inverse-OPUC analogue**,
with prescribed equispaced eigenphases and positive cyclic weights;
its `parameters.maximum_forward_node_error` is computed separately
for n<=64, rather than interpreting the *target* spectrum as the
measured floating-point eigenvalues.

`discrete_fractional_fourier(n,alpha,'spectral')` returns the **exact**
finite trace eigenphases `exp(-1i*pi*alpha*k/2)` with k=0,...,n-1
for odd n and k=0,...,n-2,n for even n. These are the eigenphases
of the Candan unitary transform (not the discrete oscillator energies);
its DFT endpoint at alpha=1 is verified using the joint parity basis.

`rational_arnoldi_circle(n,poles,nodes,weights,'spectral')` returns the
EXACT finite eigenphase nodes and uniform 1/n **trace** weights, separately
from the supplied positive e1 **cyclic** weights. The auxiliary
n-by-(n-1) Hessenberg recurrence pencil is not itself square/unitary;
only the full finite similarity closure Q'*diag(nodes)*Q is unitary.


`gauss_jacobi_transform(n,alpha,beta,'spectral')` labels the finite
orthogonal **transform** eigenphases as unknown in general.
The real quadrature nodes are eigenvalues of the Jacobi tridiagonal
recurrence, NOT of the Gauss evaluation matrix. The two FFT-fast
Chebyshev variants act on the same mathematical orthogonal matrices
as the original dense generator; no distinct spectral law is inferred.




`unitary_hessenberg_pencil` is a strict pair A,B with a
unitary quotient U=B'*A, so its generalized eigenvalues lie
on the unit circle, although their generic **trace** eigenphase
law is unknown. The projective subdiagonal ratios of A-zB
are prescribed **poles**, NOT eigenphases. The separate
`um_hessenberg_unitary_generalized` helper deforms A
without changing the pole ratios, and its spectrum is in
general **not** circle-supported.


`real_spherical_harmonic_rotation(L,alpha,beta,gamma,'spectral')`
has exactly the same finite trace phase multiset and moments
as the complex SO(3) direct sum, with exp(-i*m*delta) of
multiplicity L-|m| for all m=-(L-1),...,L-1. All eigenphases
occur in conjugate pairs as demanded by a REAL orthogonal
representation; the physical rotation angle delta is independent
of the axis and does not depend on the real/complex basis.

`gauss_hermite_transform(n,'spectral')` and
`gauss_laguerre_transform(n,alpha,'spectral')` mark the
eigenphases of the real orthogonal weighted-Gauss polynomial
evaluation matrices as UNKNOWN. Their independently available
`'quadrature'` nodes and logarithmic weights are spectral
information for the real *Jacobi tridiagonal recurrence*,
not for the unitary/orthogonal sample-evaluation matrix Q.
The additional `log_function_scaling` maps rows to physical
Hermite/Laguerre L2 function samples; it does not change
the finite spectral law of Q.


`spin_spherical_harmonic_transform(L,s,m,'spectral')`
labels the finite eigenphase law of the SQUARE fixed-(s,m)
Jacobi latitude transform as UNKNOWN. The Gauss-Jacobi
quadrature nodes are real and do not describe the
unitary/orthogonal transform eigenvalues. Its auxiliary
`'sampling'` and `'sampling_operator'` return an
L(2L-1)-by-(L²-s²) rectangular isometry rather than a
square matrix, hence there is no unitary finite trace
eigenphase counting law for that sampling matrix.
