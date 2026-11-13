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

For U of order N, the finite empirical measure is

\[
 \mu_N=\frac{1}{N}\sum_{j=1}^{N}\delta_{\arg\lambda_j},
 \qquad m_k=\frac{\operatorname{tr}(U^k)}{N}.
\]

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

The MATLAB suite tests spectral mode for all 50 public generators,
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
