function [U,properties]=cluster_floquet(n,angles,boundary,kick,mode)
%CLUSTER_FLOQUET Z2 x Z2 symmetric cluster-spin Floquet drive with pi edges.
%   U=CLUSTER_FLOQUET(N,ANGLES,BOUNDARY,KICK,MODE) has order 2^N.
%   Spin site 1 is the MOST-SIGNIFICANT tensor bit; Z|0>=|0>.
%   Two physical drive layers (up to a fixed scalar phase i^N):
%       U = P(KICK)*exp(-i*sum_j ANGLES(j)*K_j),
%       P(KICK)=prod_j X_j*exp(-i*KICK*X_j),
%       K_j=Z_(j-1)*X_j*Z_(j+1).
%   For open chains the stabilizers are j=2,...,N-1 only; for
%   periodic even chains all j=1,...,N are included (indices wrap).
%   P(KICK) is the uniform on-site pi/2+KICK X pulse, multiplied
%   by the **fixed** scalar i^N so that P(0)=prod_j X_j.
%
%   N>=4, N<=20; periodic requires even N>=4 for bipartite
%   G_odd=prod_(odd j) X_j, G_even=prod_(even j) X_j symmetries.
%   ANGLES is a real scalar or N-vector; the open endpoint angles
%   are unused (there are no endpoint cluster stabilizers).
%   Default ANGLES(j)=.43+.16*sin(.71*j), KICK=0 and BOUNDARY='open'.
%
%   At exactly KICK=0, all interior K_j commute, and in an open
%   chain the strictly local Z_1,Z_N anticommute with the Floquet U.
%   Thus every finite eigenphase lambda has an exact -lambda partner
%   (pi quasienergy separation); this is a boundary-operator statement,
%   not a claim of a CUE law. For generic nonzero KICK the on-site
%   rotation preserves both Z2 symmetries but breaks this EXACT
%   anticommutation. Finite-size splittings then need not vanish.
%   Periodic chains have no such local edge argument.
%
%   MODE='matrix' (default) yields dense unitary for N<=10.
%   MODE='operator' returns O(X,FLAG) via shared local 1/2-spin gates
%   with forward, conjugate-transpose and size actions, N<=20.
%   MODE='spectral' describes the exact finite trace eigenvalue
%   multiset and moments at the solvable KICK=0 point, otherwise
%   marks the finite analytic law unknown (not inferred using eig).
%
%   References (analogous Floquet-SPT and cluster-edge constructions;
%   the exact commuting-cluster specialization and its moment formulas
%   are derived directly, not quoted verbatim from these papers):
%   [135] A. C. Potter, T. Morimoto, A. Vishwanath, PRX 6 (2016)
%         041001. https://doi.org/10.1103/PhysRevX.6.041001
%   [136] D. V. Else, C. Nayak, PRB 93 (2016) 201103.
%         https://doi.org/10.1103/PhysRevB.93.201103
%   [137] Digital quantum simulation of Floquet symmetry-protected
%         topological phases, Nature 607 (2022) 468-473.
%         https://doi.org/10.1038/s41586-022-04854-3
properties={'square','unitary','normal','complex','parameter-dependent','scalable'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','finite','>=',4,'<=',20});
if nargin<2 || isempty(angles),angles=.43+.16*sin(.71*(1:n));end
if nargin<3 || isempty(boundary),boundary='open';end
if nargin<4 || isempty(kick),kick=0;end
if nargin<5 || isempty(mode),mode='matrix';end
validateattributes(angles,{'numeric'},{'vector','real','finite'});
if isscalar(angles),angles=repmat(angles,n,1);end
if numel(angles)~=n
    error('cluster_floquet:Angles','ANGLES must be scalar or length N.');
end
angles=angles(:);
if ~(ischar(boundary)||(isstring(boundary)&&isscalar(boundary)))
    error('cluster_floquet:Boundary','BOUNDARY must be open or periodic.');
end
boundary=lower(char(boundary));
if ~any(strcmp(boundary,{'open','periodic'}))
    error('cluster_floquet:Boundary','BOUNDARY must be open or periodic.');
end
if strcmp(boundary,'periodic') && mod(n,2)~=0
    error('cluster_floquet:Boundary','Periodic Z2 x Z2 chains require even N.');
end
validateattributes(kick,{'numeric'},{'scalar','real','finite'});
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('cluster_floquet:Mode','MODE must be text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','operator','spectral'}))
    error('cluster_floquet:Mode','Unknown output mode.');
end
if strcmp(mode,'spectral')
    U=um_cluster_floquet_spectral(n,angles,boundary,kick);
    return
end
if strcmp(mode,'matrix') && n>10
    error('cluster_floquet:ExplicitOrder','Dense output restricted to N<=10.');
end
CZ=diag([1,1,1,-1]);
gates=struct('gate',{},'sites',{});
% First CZ layer conjugates the interior X rotations into cluster K_j.
for j=1:n-1
    gates(end+1)=struct('gate',CZ,'sites',[j,j+1]); %#ok<AGROW>
end
if strcmp(boundary,'periodic')
    gates(end+1)=struct('gate',CZ,'sites',[n,1]);
    sites=1:n;
else
    sites=2:n-1;
end
for j=sites
    theta=angles(j);
    R=[cos(theta),-1i*sin(theta);-1i*sin(theta),cos(theta)];
    gates(end+1)=struct('gate',R,'sites',j); %#ok<AGROW>
end
for j=1:n-1
    gates(end+1)=struct('gate',CZ,'sites',[j,j+1]); %#ok<AGROW>
end
if strcmp(boundary,'periodic')
    gates(end+1)=struct('gate',CZ,'sites',[n,1]);
end
% X exp(-i kick X) on every site equals i exp(-i(pi/2+kick)X);
% the fixed scalar i^N makes the ideal KICK=0 pulse P=prod X.
pulse=[-1i*sin(kick),cos(kick);cos(kick),-1i*sin(kick)];
for j=1:n
    gates(end+1)=struct('gate',pulse,'sites',j); %#ok<AGROW>
end
O=@(X,flag)um_circuit_apply(X,flag,n,gates,1);
if strcmp(mode,'operator')
    U=O;
else
    U=O(eye(2^n),'notransp');
end
end
