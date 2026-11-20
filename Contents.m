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
%   unitary_rank_one_update Multiplicative unitary rank-one Hessenberg update.
%   etna_hessenberg     Type I/II/III reproducible unitary Hessenberg benchmarks.
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
%   discrete_sine        Orthonormal DST-I, II, III, IV transforms.
%   walsh_hadamard       Natural/sequency ordered fast Walsh transform.
%   dpss_transform       Explicit Slepian concentration eigenbasis.
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
%   sphere_geodesic_transport  SO(n) lift of tangent geodesic parallel transport.
%
% Orthogonal polynomial and spherical-harmonic transforms:
%   gauss_recurrence_transform  Orthogonal Gaussian transform from any positive Jacobi recurrence.
%   gauss_jacobi_transform   Jacobi/Legendre/Chebyshev/Gegenbauer/associated Legendre.
%   gauss_hermite_transform  Hermite Gauss-quadrature orthogonal basis.
%   gauss_laguerre_transform Generalized Laguerre quadrature basis.
%   wigner_rotation         Irreducible SO(3) unitary representation.
%   spherical_harmonic_rotation Direct sum of Wigner blocks.
%   periodic_sinc_shift     Real circulant fractional shift, FFT.
%
% Geometric and symplectic integration:
%   geometric_oscillator  Periodic symplectic oscillator; exact phases and metric pencil.
%   time_skew_magnus    Noncommuting time-dependent orthogonal midpoint/Cayley flow.
%   geometric_fem_wave    Finite-element Dirichlet wave; mass-weighted palindromic QEP.
%
% Nonlinear unitary dynamics:
%   lohe_matrix_flow     Coupled Lohe oscillators on U(d), geometric midpoint integration.
%
% Quantum maps, walks and algorithms:
%   split_step_walk      Periodic real-orthogonal two-state quantum walk.
%   szegedy_walk        Markov-chain product-of-reflections walk, exact spectrum.
%   coined_graph_walk   Flip-flop coins on an arbitrary weighted graph.
%   metric_quantum_graph  Bond-scattering U(k), arbitrary metric edge lengths.
%   unitary_almost_mathieu  Magnetic quasiperiodic two-state quantum walk.
%   dirac_automaton      Periodic 1D Dirac cellular automaton (size 2n).
%   grover_operator      Grover reflections with exact degeneracies.
%   quantum_baker        Balazs-Voros / Saraceno quantum baker map.
%   quantum_cat_map      Finite metaplectic quantization of hyperbolic toral map.
%   kicked_harper        Finite-torus kicked Harper Floquet map.
%   kicked_top           Spin-j quantum kicked top.
%   extended_lipkin_propagator  Nuclear two-level boson phase-transition evolution.
%   vibron_2d_esqpt      U(3) vibron ESQPT with angular-momentum sectors.
%   agassi_propagator     Fermionic nuclear pairing and particle-hole unitary dynamics.
%   floquet_ising        Integrable kicked transverse-field Ising chain.
%   selfdual_kicked_ising Self-dual kicked Ising with longitudinal fields.
%   driven_kitaev        Two-step driven Kitaev BdG operator.
%   kicked_rotor         Quantum kicked rotor (Fourier grid and kinetic phases).
%   dual_unitary_circuit Dual-unitary two-qubit brickwork Floquet circuit.
%   xxz_floquet          Integrable six-vertex XXZ brickwork Floquet circuit.
%   kicked_bose_hubbard Two-mode bosonic Floquet dimer, order n+1.
%   jaynes_cummings     Sparse cavity-QED propagator, order 2(n+1).
%   chalker_coddington  Directed sparse quantum-Hall network on a torus.
%
% Matrix functions, polar factors and unitary extensions:
%   cayley_jacobi        Cayley image of a Jacobi recurrence matrix.
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
%
% External nonlinear/palindromic benchmarks: metadata/eigenvalue_problems.json
% Public OPUC helpers (not Anymatrix generators):
%   szego_jacobi_measure     Normalized Jacobi measure and Szego lift.
%   szego_bernstein_measure  Explicit Poisson/Bernstein-Szego measure.
%   szego_circular_jacobi_measure  Generalized Jacobi circle weight.
%   szego_rogers_measure   Wrapped-normal measure.
%   szego_geronimus_measure Constant-coefficient OPUC measure.
%   szego_periodic_measure Periodic Verblunsky OPUC measure.
%   szego_bernstein_finite_measure General rational Bernstein-Szego weight.
%   szego_verblunsky        Measure -> interior Schur coefficients.
%   szego_matrix            Measure -> finite multiplication matrix.
%   szego_quadrature        e1 spectral nodes and positive weights.
%   opuc_christoffel         Positive finite-circle measure transform, inverse CMV.
%
% Geometry, geodesics, and parallel transport:
%   grassmann_geodesic        Grassmann geodesic ambient unitary lift.
%   geodesic_holonomy         Closed spherical polygon parallel transport.
%   connection_cycle          Sparse cyclic discrete-connection transport.
%   bundle_transport          General block-permutation fiber holonomies.
%   lie_geodesic_transport    Orthogonal Lie-algebra geodesic transport.
%
%   stiefel_geodesic         Canonical Stiefel geodesic ambient lift.
%   grassmann_holonomy       Non-Abelian tautological geodesic-loop holonomy.
%
%   hyperbolic_triangle_holonomy  Poincare disk triangle, SO(2) holonomy.
%   lattice_wilson_loop          Non-Abelian U(n) gauge-link products.
%
% Test runner: test_run. See docs/MATRIX_ENCYCLOPEDIA.md for references.
