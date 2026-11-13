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
%   finite_cmv           Five-diagonal CMV with specified coefficients.
%   opuc_matrix          General finite unitary twisted-core matrix.
%   jacobi_szego         Jacobi-Szego/Geronimus family, four representations.
%   bernstein_szego      Poisson/Bernstein-Szego measure, four representations.
%   circular_beta        Sparse random CMV with circular-beta eigenvalue law.
%   circular_jacobi      Generalized circular Jacobi measure (complex parameters).
%   rogers_szego         Wrapped Gaussian / Rogers-Szego orthogonal CMV.
%   geronimus_cmv        Constant real or complex Verblunsky coefficients.
%   periodic_cmv         Periodic real or complex Verblunsky coefficients.
%   bernstein_szego_finite  General finite Bernstein-Szego rational circle weight.
%
% Random and signal-processing:
%   haar_unitary         Haar random U(n) via phase-corrected complex QR.
%   haar_orthogonal      Haar random O(n) via sign-corrected real QR.
%   haar_wavelet         Orthogonal finite Haar wavelet analysis transform.
%
% Quantum maps, walks and algorithms:
%   split_step_walk      Periodic real-orthogonal two-state quantum walk.
%   dirac_automaton      Periodic 1D Dirac cellular automaton (size 2n).
%   grover_operator      Grover reflections with exact degeneracies.
%   quantum_baker        Balazs-Voros / Saraceno quantum baker map.
%   kicked_harper        Finite-torus kicked Harper Floquet map.
%   kicked_top           Spin-j quantum kicked top.
%   floquet_ising        Integrable kicked transverse-field Ising chain.
%   selfdual_kicked_ising Self-dual kicked Ising with longitudinal fields.
%   driven_kitaev        Two-step driven Kitaev BdG operator.
%   kicked_rotor         Quantum kicked rotor (Fourier grid and kinetic phases).
%   dual_unitary_circuit Dual-unitary two-qubit brickwork Floquet circuit.
%   xxz_floquet          Integrable six-vertex XXZ brickwork Floquet circuit.
%
% Matrix functions, polar factors and unitary extensions:
%   cayley_jacobi        Cayley image of a Jacobi recurrence matrix.
%   polar_unitary        Unitary polar factor of nonsingular matrix.
%   julia_unitary        Finite Julia/Halmos unitary of a contraction.
%
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
%
% Test runner: test_run. See docs/MATRIX_ENCYCLOPEDIA.md for references.
