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
%   kicked_top           Spin-j Floquet top, diagonal torsion times Jy rotation.
%   floquet_ising        Integrable kicked Ising, 2^n-dimensional spin chain.
%   selfdual_kicked_ising Self-dual disordered kicked Ising Floquet operator.
%   driven_kitaev        Two-step Kitaev chain in 2n-dimensional BdG form.
%
% Matrix functions, polar factors and unitary extensions:
%   cayley_jacobi        Cayley image of a Jacobi recurrence matrix.
%   polar_unitary        Unitary polar factor of nonsingular matrix.
%   julia_unitary        Finite Julia/Halmos unitary of a contraction.
%
% Public OPUC helpers (not Anymatrix generators):
%   szego_jacobi_measure     Normalized Jacobi measure and Szego lift.
%   szego_bernstein_measure  Explicit Poisson/Bernstein-Szego measure.
%   szego_verblunsky        Measure -> interior Schur coefficients.
%   szego_matrix            Measure -> finite multiplication matrix.
%   szego_quadrature        e1 spectral nodes and positive weights.
%
% Test runner: test_run. See docs/MATRIX_ENCYCLOPEDIA.md for references.
