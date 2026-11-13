%QUICKSTART Minimal OPUC and spectral-band examples.
% Run from the repository root after: addpath(pwd)
%
% A finite unitary matrix has an atomic eigenvalue-counting measure.
% The continuous bands below describe an infinite/large-system limit.

% 1. Measure-first Jacobi--Szego construction and cyclic-vector weights.
mu = szego_jacobi_measure(0,0);
alpha = szego_verblunsky(mu,15);
U = szego_matrix(mu,16,'cmv',1);
H = szego_matrix(mu,16,'hessenberg',1);
[nodes,weights] = szego_quadrature(U);
assert(norm(full(U'*U-speye(16)),'fro') < 1e-12)
assert(numel(alpha)==15 && numel(nodes)==16)
assert(abs(sum(weights)-1) < 1e-12)
assert(abs(trace(U)-trace(H)) < 1e-12)
% These weights are NOT the equal weights of the empirical trace measure.

% 2. A sparse Dirac automaton with exact gapped unit-circle bands.
nSites = 48;
mass = 0.45;
D = dirac_automaton(nSites,mass);       % order 2*nSites
nu = sqrt(1-mass^2);
k = (0:nSites-1)';
omega = acos(nu*cos(2*pi*k/nSites));
edge = asin(mass);                     % positive band begins at edge
thetaD = angle(eig(full(D)));          % finite eigenphases in [-pi,pi]
% Limiting bands: [edge,pi-edge] and [-pi+edge,-edge].
assert(all(abs(thetaD)>=edge-1e-12))
assert(all(abs(thetaD)<=pi-edge+1e-12))
assert(abs(trace(D)-2*sum(cos(omega)))<1e-10)

% 3. Periodic CMV: band structure refers to the INFINITE periodic operator.
% The finite unitary has discrete eigenvalues and may have boundary levels
% inside infinite-operator gaps.
pattern = [0.48;-0.31];
P = periodic_cmv(48,pattern);          % finite sparse CMV
Pop = periodic_cmv(48,pattern,'cmv',1,[],'operator');
v = exp(1i*.1*(0:47)');
assert(norm(Pop(v,'notransp')-P*v)<1e-10)

fprintf('CMV: %d nodes, cyclic weights sum %.12f.\n',numel(nodes),sum(weights));
fprintf('Dirac: two bands, gap half-width %.5f radians about 0 and pi.\n',edge);
fprintf('Periodic CMV: %d finite eigenvalues; infinite bands require a discriminant.\n',size(P,1));
