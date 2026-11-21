function P=um_hessenberg_unitary_generalized(n,poles,seed,gamma,variant)
%UM_HESSENBERG_UNITARY_GENERALIZED Nonunitary RQR-type pencil benchmark.
%   P=UM_HESSENBERG_UNITARY_GENERALIZED(N,POLES,SEED,GAMMA,VARIANT)
%   returns A-z*B with A upper Hessenberg (NOT unitary in general)
%   and B unitary upper Hessenberg. The default GAMMA=0.4 produces
%       A = A0 + GAMMA*diag((1:N)/(N+1)), B=B0,
%   with A0,B0 the strict unitary pair from
%   UNITARY_HESSENBERG_PENCIL(...,'pencil').
%   Thus the subdiagonal pole ratios A(j+1,j)/B(j+1,j) are EXACTLY
%   preserved, yet neither A nor B'*A must be unitary.
%
%   This is an explicitly NONUNITARY generalized-eigenproblem helper,
%   NOT an Anymatrix unitary generator. It represents a broader
%   Hessenberg--unitary-Hessenberg pole-swapping test class.
%   P.U=B'*A, P.A, P.B, P.poles, and classification are returned.
%   Its generalized eigenvalues match eig(P.U), but need not lie
%   on the unit circle. Finite zero and infinite pole constraints
%   imply reducible factors as in the strict base construction.
%
%   [163] Camps et al., SIAM Spotlights (2025):
%   https://doi.org/10.1137/1.9781611978377
%   [164] Camps et al., ETNA 52 (2020):
%   https://doi.org/10.1553/etna_vol52s480
if nargin<1||isempty(n),n=8;end
if nargin<2,poles=[];end
if nargin<3,seed=[];end
if nargin<4||isempty(gamma),gamma=.4;end
if nargin<5||isempty(variant),variant='complex';end
validateattributes(gamma,{'numeric'},{'scalar','real','finite','positive'});
strict=unitary_hessenberg_pencil(n,poles,seed,variant,'pencil');
B=strict.B;
A=strict.A+gamma*diag((1:n)'/(n+1));
P=struct('A',A,'B',B,'U',B'*A,'poles',strict.poles, ...
    'target_poles',strict.target_poles, ...
    'proper',strict.proper,'gamma',gamma, ...
    'classification','generalized_hessenberg_unitary_hessenberg', ...
    'unitary_factors',false,'strict_unitary_pencil',false, ...
    'notes',['A is upper Hessenberg and generally nonunitary; ', ...
        'B is a unitary upper Hessenberg matrix. The quotient ', ...
        'B''*A is generally nonunitary, and its generalized ', ...
        'eigenvalues need not lie on the circle. This helper ', ...
        'MUST NOT be registered as an Anymatrix unitary family.']);
end
