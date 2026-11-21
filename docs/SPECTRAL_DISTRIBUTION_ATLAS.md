# Spectral Distribution Atlas for Unitary Matrices

**Status:** Classification draft (21 November 2026). This is a taxonomy and a set of candidates for a future illustrated atlas. No figures are claimed to be present yet. The classification distinguishes proven/implemented analytic references from proposed or exploratory cases.

## Purpose and mathematical conventions

The atlas groups the spectral **measures and density shapes** occurring among the generators of [`anymatrix-unitary`](../README.md), rather than grouping them by application. The same matrix can belong to several morphological classes. For instance, the large-size Dirac automaton trace law is simultaneously absolutely continuous, two-band, gapped, conjugation symmetric, and singular at band edges.

For a unitary $N\times N$ matrix $U_N$, write $\lambda_j=e^{i\theta_j}$, $\theta_j\in[-\pi,\pi)$. The finite *normalized trace counting measure* is

$$
\mu_N=\frac1N\sum_{j=1}^N\delta_{\theta_j},\qquad
m_k=\int e^{ik\theta}\,d\mu_N(\theta)=\frac{\operatorname{tr}(U_N^k)}{N}.
$$

**Every finite-dimensional measure $\mu_N$ is atomic.** Continuous distributions in this document refer to (i) a verified or derived large-$N$ weak limit of $\mu_N$, (ii) an ensemble average, or (iii) an infinite-dimensional operator's spectral measure at a distinguished vector. These are different quantities and **must not** be exchanged without justification. In particular, a finite CMV matrix's uniform trace weights differ from its Gauss–Szegő quadrature weights and from the infinite OPUC measure.

For a general Borel probability measure $\mu$ on the circle, the Lebesgue decomposition is

$$
 d\mu(\theta)=\rho_{\rm ac}(\theta)\,d\theta+\sum_j w_j\,\delta_{\theta_j}+d\mu_{\rm sc}(\theta).
$$

Here $\mu_{\rm sc}$ is **singular continuous**: it has no atoms but is singular relative to angular Lebesgue measure. The nonnegative density $\rho_{\rm ac}$ need not integrate to one when atoms or a singular-continuous component are present. Throughout, a density is defined **with respect to $d\theta$**, not $d\theta/(2\pi)$; the uniform density is $1/(2\pi)$.

The repository's authoritative interfaces are [`docs/SPECTRAL_REFERENCE.md`](SPECTRAL_REFERENCE.md), [`unitary_spectral_reference.m`](../unitary_spectral_reference.m), and [`metadata/catalog.json`](../metadata/catalog.json). Descriptors `S.finite`, `S.limit`, `S.ensemble`, and `S.underlying` must remain distinct.

## Classification system: several independent axes

### Axis 1 — Type and provenance of the measure

| Tag | Definition | Notes |
|---|---|---|
| `finite:atomic` | Normalized trace measure of any specified finite unitary matrix. | Universal at finite $N$; multiplicities matter. |
| `limit:ac` | Large-$N$ trace law with density $\rho(\theta)$. | Requires a defined sequence of matrices and a justified limit. |
| `limit:atomic` | Large-$N$ trace law with one or more surviving atoms. | High multiplicities can persist after normalization. |
| `limit:mixed` | Limit with both atomic and continuous parts. | Can also have singular-continuous mass. |
| `ensemble:ac` | Expected one-point law of a random ensemble. | Not the density of a finite realization. |
| `underlying:ac` / `underlying:mixed` | Infinite-dimensional cyclic-vector (OPUC/CMV) spectral law. | Not automatically a limiting trace law of finite truncations. |
| `singular-continuous` | Continuous measure singular to Lebesgue measure. | Needs a rigorous spectral-type model; a finite plot cannot establish it. |
| `unknown` | No analytic reference or rigorously identified spectral type. | Do not infer one solely from a finite sample. |

### Axis 2 — Geometry of spectral support

| Tag | Description | What an atlas image should reveal |
|---|---|---|
| `finite-points` | Support at finitely many positions on $\mathbb T$. | Distinct phases, multiplicities, and separation. |
| `full-circle` | Support closure is the entire unit circle. | Uniformity or nonuniformity of angular mass. |
| `one-band` | Support is one proper closed arc. | Two endpoints and the complementary gap. |
| `multi-band` | Support is two or more proper arcs. | Band endpoints and the number/size of genuine gaps. |
| `isolated-atom` | An atom lies away from the continuous support. | A discrete phase inside an open gap. |
| `thin-or-Cantor` | Support has a Cantor-like/fractal structure. | Hierarchical gaps under scale refinement; no claim of density from a finite histogram. |

On the circle, arcs meeting across $-\pi$ and $\pi$ are **connected**; the choice of angle cut must not create an artificial second band. Likewise, an empty region between a finite sample's eigenvalues is not automatically a spectral gap of its limiting support.

### Axis 3 — Local density behavior

| Tag | Analytic feature | Visualization concern |
|---|---|---|
| `constant` | $\rho\equiv 1/(2\pi)$ on the full circle. | Flat baseline. |
| `smooth-peaked` | Smooth positive, nonuniform density; possibly multimodal. | Peak width and small secondary modes. |
| `hard-gap` | Exactly zero density on an open arc. | Smoothing fills the gap artificially. |
| `inverse-sqrt-edge` | $\rho(\theta)\sim C/\sqrt{|\theta-\theta_*|}$ near a band edge, on the occupied side. | Integrable divergence; evaluate analytic reference away from endpoints. |
| `soft-edge` | $\rho\to0$ toward a band endpoint, with some positive power. | Peak-free end of support, unlike inverse-square-root. |
| `cusp-or-power-law` | A zero, divergence, or nondifferentiability at an interior/special phase. | Distinguish algebraic behavior from a narrow smooth peak. |
| `multimodal` | Two or more separated local maxima, with or without gaps. | Multiple peaks need not mean separate support bands. |
| `localized-cluster` | Finite samples concentrate strongly in small angular regions. | Finite-size feature; do not assert an atom in the limit. |
| `degenerate` | Exactly repeated phases. | Count or annotate multiplicity; coincident hairlines are invisible. |

A genuine point mass is **not** an infinite value of $\rho_{\rm ac}$; it belongs in the atomic component. Conversely, a narrow but smooth Poisson/Bernstein–Szegő peak is not an atom.

### Axis 4 — Additional spectral organization

- `conjugate-symmetric`: eigenphases occur at $\pm\theta$ for real orthogonal matrices, with possible eigenvalues $\pm1$.
- `pi-paired`: special drive/symmetry models can pair phases separated by $\pi$; verify for the specific model and sector.
- `wrapped` / `aliased`: $\theta$ arises by reduction of a real dispersion or quasienergy modulo $2\pi$.
- `sector-dependent`: pooled and individual symmetry sectors can have different multiplicities and apparent density.
- `parameter-transition`: band closure, gap opening, an atom entering/leaving the spectrum, or a change of spectral type under parameter variation.
- `level-correlated`: a statistical property **not determined by the one-point density** (e.g. CUE repulsion versus independently distributed uniform points).

## Curated gallery candidates already represented in the repository

The entries below identify **which measure is being classified** and how strong the evidence is. These are not claims that every candidate has a closed-form limiting trace density.

| Candidate / concrete choice | Measure being classified | Tags | Evidence and caution |
|---|---|---|---|
| [`discrete_fourier`](../discrete_fourier.m), $N=64$ | Finite trace | `finite:atomic`, `finite-points`, `degenerate`, four locations $\{1,-1,\pm i\}$ | Exact reference in `unitary_spectral_reference.m`. Never label its counting law uniform. |
| [`grover_operator`](../grover_operator.m), $N=128$, marked $M=2$ | Finite trace | `finite:atomic`, `degenerate`, large dominant atom, two rare phases | Exact $\pm2\arcsin\sqrt{M/N}$ and $\pm1$ multiplicities. Useful to test loss of small masses under smoothing. |
| [`diagonal_phase`](../diagonal_phase.m) with explicitly prescribed angular clusters | Finite trace (synthetic benchmark) | `localized-cluster`, `multimodal`, optionally `hard-gap` **in a chosen limit** | Exact prescribed finite phases; no continuum limit unless a sequence is defined. |
| [`cycle_shift`](../cycle_shift.m), $N\to\infty$ | Limiting trace | `limit:ac`, `full-circle`, `constant` | Exact $N$th roots of unity; weak limit uniform. |
| [`haar_unitary`](../haar_unitary.m) / [`circular_beta`](../circular_beta.m) | Ensemble-mean density | `ensemble:ac`, `full-circle`, `constant`, `level-correlated` | Uniform expectation; individual finite samples remain atomic and have correlations. |
| [`haar_special_unitary`](../haar_special_unitary.m), finite $N$ | Ensemble-mean density | `ensemble:ac`, `full-circle`, finite-size angular modulation | Repository cites Nishigaki's exact finite-$N$ expectation; not the same law as unrestricted Haar U($N$). |
| [`bernstein_szego`](../bernstein_szego.m), $a=0.6$ | Underlying infinite OPUC vector measure | `underlying:ac`, `full-circle`, `smooth-peaked` | Poisson-type positive rational circle density. **Not** generally the finite-CMV trace density. |
| [`rogers_szego`](../rogers_szego.m), $q=0.6$ | Underlying infinite OPUC vector measure | `underlying:ac`, `full-circle`, `smooth-peaked` | Wrapped Gaussian-type law; not a finite trace limit claim. |
| [`cycle_propagator`](../cycle_propagator.m), $t=0.8$ | Limiting trace | `limit:ac`, `one-band`, `hard-gap`, `inverse-sqrt-edge` | Wrapped arcsine law; for $2|t|<\pi$ its occupied arc is $[-2|t|,2|t|]$. |
| [`cycle_propagator`](../cycle_propagator.m), $t>\pi/2$ | Limiting trace | `limit:ac`, `wrapped`, `multimodal`; potentially overlapping branch singularities | The same arcsine distribution folds modulo $2\pi$; support and edge locations depend on $t$. |
| [`dirac_automaton`](../dirac_automaton.m), $0<m<1$, e.g. $m=0.55$ | Limiting trace | `limit:ac`, `multi-band` (two arcs), `hard-gap`, `inverse-sqrt-edge`, `conjugate-symmetric` | Exact two-band density derived from the dispersion relation; at $m=0$ uniform and at $m=1$ two atoms. |
| [`spherical_harmonic_rotation`](../spherical_harmonic_rotation.m), rotation angle $\delta=c/L$, $0<c<\pi$ | Limiting trace for **size-dependent** parameters | `limit:ac`, `one-band`, `soft-edge`, piecewise linear | Triangular $\rho(\theta)=(1-|\theta|/c)/c$ on $|\theta|\le c$. Not a limit claim for fixed Euler angles. |
| [`jacobi_szego`](../jacobi_szego.m), vary Jacobi exponents $a,b>-1$ | Underlying infinite OPUC vector measure | `underlying:ac`, `full-circle`, `soft-edge` / `cusp-or-power-law` at $0,\pi$ | Explicit symmetric Szegő lift. For exponent choices, density can vanish or diverge integrably at special phases. |
| [`circular_jacobi`](../circular_jacobi.m), $a\ne0$ | Underlying infinite OPUC vector measure | `underlying:ac`, `full-circle`, `cusp-or-power-law`; asymmetric for nonzero $\eta$ | Generalized circular Jacobi weight. A density zero at one phase does **not** imply an open gap. |
| [`geronimus_cmv`](../geronimus_cmv.m), constant coefficient $a=0.3+0.2i$ | Underlying infinite OPUC vector measure | `underlying:mixed`, `one-band`, `hard-gap`, `isolated-atom`, `soft-edge` (depending on location/parameters) | `szego_geronimus_measure.m` records an arc density and an atom if $\Re a+|a|^2>0$. **Not** a finite trace-density assertion. |
| [`periodic_cmv`](../periodic_cmv.m), nonconstant period | Infinite periodic background / support | `multi-band`, `hard-gap`, parameter transitions | Known periodic spectral-band support; generic finite trace density is **not** supplied analytically by the repository. Finite CMV boundary conditions matter. |
| [`orthogonal_blocks`](../orthogonal_blocks.m), custom angles | Finite trace | `conjugate-symmetric`, optional `degenerate` or `localized-cluster` | Exact paired phases and possibly $\pm1$. Symmetry is an additional tag, not a measure type. |
| [`unitary_almost_mathieu`](../unitary_almost_mathieu.m) | Finite trace; infinite model needs separate analysis | `unknown` at finite generic parameters; research target for fractal/singular measures | Finite periodic-window closure does **not** automatically inherit infinite irrational-flux spectral-type statements. |

### Analytic shapes worth illustrating

1. **Uniform:** $\rho(\theta)=1/(2\pi)$.
2. **Smooth concentrated Poisson weight:** $\rho(\theta)=(1-|a|^2)/(2\pi|e^{i\theta}-\bar a|^2)$, $|a|<1$ (a spectral law, not the plot's external smoothing kernel).
3. **Single-arc inverse-square-root:** for the cycle propagator with $a=2|t|<\pi$, $\rho(\theta)=\bigl[\pi\sqrt{a^2-\theta^2}\bigr]^{-1}$ for $|\theta|<a$, zero elsewhere.
4. **Two-arc inverse-square-root:** for the Dirac automaton with $\nu=\sqrt{1-m^2}$, $\rho(\theta)=|\sin\theta|/[2\pi\sqrt{\nu^2-\cos^2\theta}]$ when $|\cos\theta|<\nu$, zero elsewhere.
5. **One-arc soft edges:** the special triangular trace limit $\rho(\theta)=(1-|\theta|/c)_+/c$.
6. **Algebraic zero or divergence:** the circular Jacobi underlying weight with $\eta=0$ has $\rho(\theta)\propto|1-e^{i\theta}|^{2a}$, $a>-1/2$.
7. **Mixed continuous plus atom:** the constant-Verblunsky underlying Geronimus measure, with an explicit atom when its parameter satisfies the condition above.

## Genuine gaps not yet covered by a verified database reference

### Singular continuous and hierarchical gap structure

Infinite-dimensional CMV and quantum-walk literature contains examples of purely singular-continuous measures, including limit-periodic constructions. These should form an eventual atlas section, but **none of the current finite-matrix generators should be assigned such a limiting type merely because its plot looks fractal**. We need (i) a precisely specified infinite model with theorem-backed spectral type and (ii) faithful structure-preserving finite approximants, explicitly separated from the limiting object. An interesting secondary distinction is *Cantor support versus singular-continuous measure*: neither implies the other in full generality.

### Additional deliberately synthesized stress tests

Use the existing `diagonal_phase` or `phase_circulant` generators to construct, with documented size-dependent sampling, (i) several sharply separated clusters, (ii) an absolutely continuous background plus a macroscopic atom, (iii) near-colliding atoms, (iv) uniform angular density on one or multiple arcs, and (v) continuous density plus a small isolated mass. These are valuable reference shapes but must be marked **synthetic**, not as established application-derived limit laws.

## Figure specification for the eventual illustrated atlas

Our current experimental style is the standard default, unless modified for a particular mathematical feature:

- **Axes:** complex plane with a light polar grid extending beyond the unit circle and the outer density curve; centered radial eigenphase **hairlines** straddling $|z|=1$.
- **Density:** a separate outward contour, $R(\theta)=1+s\,\rho_r(\theta)$, with its radial displacement labeled as a *visualization scale* rather than eigenvalue modulus.
- **Kernel:** periodic Poisson smoothing with default **$r=0.94$**, $0<r<1$, using $P_r(\phi)=(1-r^2)/[2\pi(1-2r\cos\phi+r^2)]$; the parameter and bandwidth must be displayed in captions.
- **Palette:** monochrome/gray, high-resolution PNG plus vector SVG or PDF; avoid dependence on color to distinguish spectra and density curves.
- **Discrete atoms:** annotate exact multiplicity/weight; several coincident eigenphases **cannot** be represented faithfully by a single unweighted hairline.
- **Reference data:** use the built-in `'spectral'` interface when available, and record the specific reference type (`finite`, `limit`, `ensemble`, or `underlying`).
- **Comparison discipline:** preserve the same radial scale when comparing smoothing parameters for the **same** example; do not interpret independently normalized radial panels as amplitude-comparable.
- **Singular continuous:** do not draw a conventional density curve as if one existed. Plot phase counts at several resolutions, cumulative distribution, or scale-dependent Poisson-smoothed approximants explicitly labeled as approximations.

### Reproducible spectral-reference starting points

```matlab
Sdft   = discrete_fourier(64, 'spectral');
Sdirac = dirac_automaton(256, 0.55, 'spectral');
Scycle = cycle_propagator(256, 0.8, 'spectral');
Shaar  = haar_unitary(144, [], 'spectral');

z_dft   = Sdft.finite.eigenvalues();
z_dirac = Sdirac.finite.eigenvalues();
z_cycle = Scycle.finite.eigenvalues();
% Shaar.ensemble.density is the exact *mean*, not a sampled spectrum.
```

The eventual renderer should take these reference descriptors as inputs; for a Haar sample, construct a seeded `haar_unitary(144, seed)` and use `eig` on that particular sample if the finite eigenphases are desired.

The atlas should first offer a **one-page classification plate** with representative examples, followed by one section per class with parameter values, a reproducible generator invocation, analytic reference formulas where known, and cautions about spectral type. Actual rendered figures are deferred to the next stage of this draft.

## Further primary references

1. B. Simon, *Orthogonal Polynomials on the Unit Circle*, Parts 1–2, AMS, 2005. [Part 1](https://doi.org/10.1090/coll054.1), [Part 2](https://doi.org/10.1090/coll/054.2).
2. M. J. Colbrook, “Computing spectral measures and spectral types,” *Commun. Math. Phys.* **384** (2021), 433–501. [DOI](https://doi.org/10.1007/s00220-021-04072-4). Includes continuous-and-atomic Geronimus illustrations.
3. J. Fillman and D. C. Ong, “Purely singular continuous spectrum for limit-periodic CMV operators with applications to quantum walks,” *J. Funct. Anal.* **272** (2017), 5107–5143. [Preprint](https://arxiv.org/abs/1610.06159).
4. A. Bisio, G. M. D’Ariano, and A. Tosini, “Quantum field as a quantum cellular automaton: The Dirac free evolution in one dimension,” *Ann. Phys.* **354** (2015), 244–264. [DOI](https://doi.org/10.1016/j.aop.2014.12.016).
5. F. Mezzadri, “How to generate random matrices from the classical compact groups,” *Notices AMS* **54** (2007), 592–604. [arXiv](https://arxiv.org/abs/math-ph/0609050).
6. S. Nishigaki, *Prog. Theor. Exp. Phys.* (2024), 021A02, for nonuniform exact finite-$N$ SU and unimodular circular ensemble mean densities. [DOI](https://doi.org/10.1093/ptep/ptae018).

**Editorial rule:** never promote an exact *finite* spectrum, a reference for an *underlying OPUC measure*, or an *ensemble mean* into a claim of a *limiting normalized-trace density* without an additional theorem or derivation.