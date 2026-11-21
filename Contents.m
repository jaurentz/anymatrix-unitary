% UNITARY  Unitary and orthogonal test matrix generators (Anymatrix group).
%
% Prescribed spectra / circulants:
%   cycle_shift          Sparse cyclic permutation with exact spectrum.
%   weighted_shift       Sparse phase-weighted cyclic permutation.
%   diagonal_phase       Sparse unitary with prescribed eigenphases.
%   phase_circulant      Fourier-conjugated prescribed phase matrix.
%   cycle_propagator     Exactly solved tight-binding propagator on a cycle.
%   orthogonal_blocks    Prescribed real orthogonal rotation blocks.
%
% Orthogonal polynomials on the unit circle:
%   finite_cmv           Five-diagonal CMV; terminal 'pencil' mode for factors A,B.
%   inverse_opuc         Inverse eigenphases/positive cyclic weights to finite CMV.
%   schur_hessenberg_flow Isospectral Hessenberg QR/Schur flow, evolving cyclic weights.
%   discrete_schur_cmv   Positive finite CMV weight-step spectral trajectory.
%   abl_ladik_cmv       Exact finite defocusing Ablowitz-Ladik Lax isospectral CMV.
%   circle_polynomial_unitary Unit-circle polynomial to strict unitary core realization.
%   opuc_matrix          General finite unitary twisted-core matrix.
%   jacobi_szego         Jacobi-Szego/Geronimus family, four representations.
%   bernstein_szego      Poisson/Bernstein-Szego measure, four representations.
%   circular_beta        Sparse random CMV with circular-beta eigenvalue law.
%   circular_jacobi      Generalized circular Jacobi measure (complex parameters).
%   rogers_szego         Wrapped Gaussian / Rogers-Szego orthogonal CMV.
%   geronimus_cmv        Constant real or complex Verblunsky coefficients.
%   periodic_cmv         Periodic real or complex Verblunsky coefficients.
%   bernstein_szego_finite  General finite Bernstein-Szego rational circle weight.
%   rational_opuc       Pole-parameterized finite rational OPUC unitary model.
%   rational_arnoldi_circle Full weighted rational Arnoldi basis, rectangular pencil and unitary closure.
%   unitary_rank_one_update Multiplicative unitary rank-one Hessenberg update.
%   etna_hessenberg     Type I/II/III and ETNA 2015 Examples 7.1/7.5 Hessenberg benchmarks.
%   lee_yang_ising      Sparse twisted Floquet-CMV for Ising Lee-Yang zeros.
%   scattering_zipper   General multichannel scattering zipper with boundary mixing.
%
% Random and signal-processing:
%   haar_unitary         Haar random U(n) via phase-corrected complex QR.
%   haar_orthogonal      Haar random O(n) via sign-corrected real QR.
%   haar_special_unitary Haar random determinant-one SU(n) ensemble.
%   haar_unimodular_coe Haar symmetric SU(n)/SO(n) circular ensemble.
%   haar_unimodular_cse Haar self-dual SU(2n)/Sp(2n) circular ensemble.
%   haar_tensor_power    Correlated repeated Haar tensor factors.
%   fermionic_exterior_power  Number-conserving exterior-power unitary lift.
%   bosonic_symmetric_power  Bosonic Fock-sector symmetric-power lift.
%   optical_mesh            Lossless triangular/rectangular photonic coupler meshes.
%   haar_block_tensor_power Two independent Haar factors in tensor blocks.
%   tensor_permutation  Orthogonal tensor-factor permutation matrix.
%   entangled_reflection Maximally entangled rank-one unitary phase reflector.
%   haar_wavelet         Orthogonal finite Haar wavelet analysis transform.
%   daubechies_wavelet   Periodized db2/db3/db4 or custom QMF; wavelet-packet trees.
%
% Non-quantum signal-processing transforms:
%   orthogonal_dct       Orthonormal DCT-I, II, III, IV finite transforms.
%   discrete_hartley     Symmetric real orthogonal Hartley involution.
%   cazac_circulant      Zadoff-Chu circulant complex Hadamard matrices.
%   fresnel_transform    Discrete Fresnel chirp transform (OCDM).
%   affine_fourier       AFDM chirp-DFT-chirp modulation transform.
%   finite_zak           Finite Zak transform on a factored signal length.
%   time_frequency_shift Sparse Weyl-Heisenberg modulation/translation.
%   periodic_hilbert     Real orthogonal periodic Hilbert completion.
%   discrete_fourier     Normalized DFT and exact fourth-root spectrum.
%   discrete_fractional_fourier Candan DFT-compatible fractional Fourier (Hermite parity basis).
%   discrete_sine        Orthonormal DST-I, II, III, IV transforms.
%   walsh_hadamard       Natural/sequency ordered fast Walsh transform.
%   dpss_transform       Explicit 1D Slepian concentration eigenbasis.
%   spherical_slepian    Real spherical polar-cap Slepian orthogonal eigenbasis.
%   spherical_lune_slepian Exact non-axisymmetric meridional-lune Slepian eigenbasis.
%   graph_fourier        Explicit weighted-graph Laplacian eigenbasis.
%   blaschke_potapov    Rational inner / lossless rank-one product frequency response.
%   paraunitary_fir     Multichannel FIR polyphase, exact finite periodic unitary.
%   orthogonal_butterfly Recursive real orthogonal butterfly (fast implicit action).
%   complex_hadamard     Dita-deformed Fourier-block complex Hadamard matrices.
%   butson_hadamard      Integer-phase Tao BH(6,3) and Paley I BH(p+1,2).
%   krawtchouk_transform Orthogonal binomial Krawtchouk polynomial transform.
%   hahn_transform        Orthogonal Hahn transform via Jacobi difference operator.
%   dual_hahn_transform   Dual Hahn quadratic-grid polynomial transform.
%
% Riemannian geometry:
%   hyperbolic_tangent_holonomy Genuine H^n tangent SO(n) Levi-Civita polygon holonomy.
%   sphere_geodesic_transport  SO(n) lift of tangent geodesic parallel transport.
%   berry_wilson_loop     Polar-overlap Berry/Wilczek-Zee from moving frames.
%   graph_connection_propagator  exp(-itL) of a general graph connection Laplacian.
%
% Orthogonal polynomial and spherical-harmonic transforms:
%   gauss_recurrence_transform  Orthogonal Gaussian transform from any positive Jacobi recurrence.
%   gauss_jacobi_transform   Jacobi transforms; true FFT operator for Chebyshev I/II.
%   gauss_hermite_transform  Hermite Gauss orthogonal basis, stable Sturm signs/log weights.
%   hermite_oscillator    Finite GH harmonic oscillator; exact fractional phase flow.
%   gauss_laguerre_transform Generalized Laguerre basis, stable Sturm signs/log weights.
%   fourier_laguerre_ball  3D radial-zonal Fourier-Laguerre, exact square tensor and rectangular sampling.
%   wigner_rotation         Irreducible SO(3) unitary representation.
%   spherical_harmonic_rotation Complex Wigner direct sum; O(L^3) cached block operator.
%   so3_peter_weyl_transform Fixed-pair orthogonal Wigner-Jacobi; rectangular SO3 FFT isometry.
%   spin_spherical_harmonic_transform Fixed spin/m orthogonal Gauss-Jacobi, rectangular FFT sampling.
%   real_spherical_harmonic_rotation Real orthogonal Condon-Shortley Wigner direct sum, block operator.
%   periodic_sinc_shift     Real circulant fractional shift, FFT.
%
% Geometric and symplectic integration:
%   geometric_oscillator  Periodic symplectic oscillator; exact phases and metric pencil.
%   gyroscopic_oscillator Conservative gyroscopic QEP, Gauss step, circle pencil.
%   moser_veselov      Implicit Moser-Veselov SO(n) rigid-body increment and trajectory.
%   time_skew_magnus    Noncommuting time-dependent orthogonal midpoint/Cayley flow.
%   geometric_fem_wave    Finite-element Dirichlet wave; mass-weighted palindromic QEP.
%   nonuniform_fem_wave   Heterogeneous/noncommuting consistent FE; metric Gauss step.
%   tensor_fem_wave       Separable 2D/3D consistent FE wave; exact phases, DST actions.
%
% Nonlinear unitary dynamics:
%   lohe_matrix_flow     Coupled Lohe oscillators on U(d), geometric midpoint integration.
%
% Quantum maps, walks and algorithms:
%   split_step_walk      Periodic real-orthogonal two-state quantum walk.
%   szegedy_walk        Markov-chain product-of-reflections walk, exact spectrum.
%   coined_graph_walk   Flip-flop coins on an arbitrary weighted graph.
%   metric_quantum_graph  Bond-scattering U(k), arbitrary metric edge lengths.
%   resonant_scattering    Open-channel flux-normalized unitary S(E), cavity poles/delay.
%   electric_quantum_walk  Finite electric two-state ring, exact cyclic CMV factors.
%   unitary_almost_mathieu  Magnetic quasiperiodic two-state quantum walk.
%   dirac_automaton      Periodic 1D Dirac cellular automaton (size 2n).
%   flat_torus_spin_dirac Exact 2D flat spin-Dirac FFT propagator and metric model.
%   grover_operator      Grover reflections with exact degeneracies.
%   quantum_baker        Balazs-Voros / Saraceno quantum baker map.
%   koopman_mpedmd        Weighted snapshot data -> Gram-Procrustes unitary Q.
%   finite_koopman_torus   Classical finite torus Koopman permutation; exact cycles.
%   quantum_cat_map      Finite metaplectic quantization of hyperbolic toral map.
%   kicked_harper        Finite-torus kicked Harper Floquet map.
%   kicked_top           Spin-j quantum kicked top.
%   extended_lipkin_propagator  Nuclear two-level boson phase-transition evolution.
%   vibron_2d_esqpt      U(3) vibron ESQPT with angular-momentum sectors.
%   agassi_propagator     Fermionic nuclear pairing and particle-hole unitary dynamics.
%   floquet_ising        Integrable kicked transverse-field Ising chain.
%   selfdual_kicked_ising Self-dual kicked Ising with longitudinal fields.
%   cluster_floquet      Z2xZ2 cluster-spin drive with exact pi edge pairing.
%   driven_kitaev        Two-step driven Kitaev BdG operator.
%   kicked_rotor         Quantum kicked rotor (Fourier grid and kinetic phases).
%   dual_unitary_circuit Dual-unitary two-qubit brickwork Floquet circuit.
%   clifford_circuit     Ordered stabilizer Clifford gates; exact atomic presets.
%   fibonacci_braid      Fibonacci anyon fusion sectors and F/R braid words.
%   xxz_floquet          Integrable six-vertex XXZ brickwork Floquet circuit.
%   kicked_bose_hubbard Two-mode bosonic Floquet dimer, order n+1.
%   jaynes_cummings     Sparse cavity-QED propagator, order 2(n+1).
%   chalker_coddington  Directed sparse quantum-Hall network on a torus.
%
% Matrix functions, polar factors and unitary extensions:
%   cayley_jacobi        Cayley image of a Jacobi recurrence matrix.
%   toda_cayley          Finite Jacobi Toda QR flow and symmetric unitary Cayley image.
%   polar_unitary        Unitary polar factor of nonsingular matrix.
%   julia_unitary        Finite Julia/Halmos unitary of a contraction.
%
% Unitary eigenvalue problems (strict generalized pencils):
%   mobius_unitary_pencil       Mobius-transformed cyclic unitary pencil.
%   unitary_delay_pencil        Lossless integer-delay secular polynomial lift.
%   block_cmv_pencil           True unitary matrix with matrix-valued CMV cores, operator and pencil.
%   floquet_factor_pencil      Two-layer split-step orthogonal pencil.
%   cayley_hermitian_pencil    Sparse Cayley pencil of a real Jacobi matrix.
%   left_weighted_unitary_pencil  Ill-scaled coefficient pair with unitary quotient.
%   unitary_hessenberg_pencil  Core-factored strict unitary Hessenberg pair with prescribed projective poles.
%
% External nonlinear/palindromic benchmarks: metadata/eigenvalue_problems.json
% Public OPUC helpers (not Anymatrix generators):
%   szego_jacobi_measure     Normalized Jacobi measure and Szego lift.
%   szego_bernstein_measure  Explicit Poisson/Bernstein-Szego measure.
%   szego_circular_jacobi_measure  Generalized Jacobi circle weight.
%   szego_rogers_measure   Wrapped-normal measure.
%   szego_geronimus_measure Constant-coefficient OPUC mixed arc/atom measure.
%   js_geronimus_density   Exact infinite e1 continuous arc density.
%   js_geronimus_poisson   Exact Poisson smoothing including an atom.
%   szego_periodic_measure Periodic Verblunsky OPUC measure.
%   szego_bernstein_finite_measure General rational Bernstein-Szego weight.
%   szego_verblunsky        Measure -> interior Schur coefficients.
%   szego_matrix            Measure -> finite multiplication matrix.
%   szego_quadrature        e1 spectral nodes and positive weights.
%   opuc_christoffel         Positive finite-circle measure transform, inverse CMV.
%
% Geometry, geodesics, and parallel transport:
%   grassmann_tangent_holonomy Intrinsic real/complex Grassmann Levi-Civita tangent transport.
%   grassmann_geodesic        Grassmann geodesic ambient unitary lift.
%   geodesic_holonomy         Closed spherical polygon parallel transport.
%   connection_cycle          Sparse cyclic discrete-connection transport.
%   bundle_transport          General block-permutation fiber holonomies.
%   lie_geodesic_transport    Orthogonal Lie-algebra geodesic transport.
%
%   stiefel_parallel_transport Intrinsic real Stiefel tangent parallel isometry.
%   real_flag_parallel_transport Normal-metric real flag tangent parallel isometry.
%   stiefel_geodesic         Canonical Stiefel geodesic ambient lift.
%   grassmann_holonomy       Non-Abelian tautological geodesic-loop holonomy.
%
%   hyperbolic_triangle_holonomy  Poincare disk triangle, SO(2) holonomy.
%   lattice_wilson_loop          Non-Abelian U(n) gauge-link products.
%
% Test runner: test_run. See docs/MATRIX_ENCYCLOPEDIA.md for references.

% Additional generalized pencil helper (NOT a unitary generator):
%   um_hessenberg_unitary_generalized  General nonunitary Hessenberg--unitary-Hessenberg pair.
