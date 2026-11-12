%SHOWCASE_NEW_FAMILIES Short examples to explore different spectral regimes.
% Run from the repository root after addpath(pwd).

% Exact prescribed eigenphase distribution with dense circulant mixing.
U=phase_circulant(32,linspace(-pi/2,pi/2,32));
assert(norm(U'*U-eye(32),'fro')<1e-11)

% Random sparse circular-beta matrix: eigenvalue law is CUE for beta=2.
C=circular_beta(64,2,42);
assert(norm(full(C'*C-speye(64)),'fro')<1e-11)

% Analytic gapped Dirac band structure and point spectrum at order 2*n.
D=dirac_automaton(40,.45);
assert(norm(full(D'*D-speye(80)),'fro')<1e-11)

% A unitary extension of a strict contraction.
T=.3*full(cycle_shift(8));
J=julia_unitary(8,T);
assert(norm(J'*J-eye(16),'fro')<1e-11)

% Poisson/Szego continuous weight -> coefficients -> CMV matrix.
mu=szego_bernstein_measure(.3+.2i);
alpha=szego_verblunsky(mu,19);
B=szego_matrix(mu,20,'cmv',1);
[z,w]=szego_quadrature(B);
assert(abs(sum(w.*z)-conj(mu.alpha0))<1e-10)

fprintf('Explored %d eigenphases in total.\n',numel(eig(U))+numel(eig(full(C)))+numel(eig(full(D)))+numel(eig(J))+numel(eig(full(B))));
