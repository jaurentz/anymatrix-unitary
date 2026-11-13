function [U,properties] = xxz_floquet(n,u,gamma,boundary,output)
%XXZ_FLOQUET Brickwork Floquet circuit of unitary trigonometric XXZ gates.
%   Uses the normalized checked six-vertex R-matrix gate G(u,gamma):
%       G=diag(a,[...],a), a=sinh(u+i*gamma)/rho,
%       middle block [i*sin(gamma),sinh(u);sinh(u),i*sin(gamma)]/rho,
%       rho=sqrt(sinh(u)^2+sin(gamma)^2), Delta=cos(gamma).
%   Full U=(even-bond gates)*(odd-bond gates), order 2^n.
%   For periodic boundaries, n even>=4 and wrap bond (n,1).
%   Default boundary: periodic for even n>=4, open otherwise.
%   For open boundaries, n>=2 (open finite circuit extension;
%   integrable boundary conditions are not claimed automatically).
%   Each gate preserves Z-magnetization; finite U commutes with
%   total Z. Yang-Baxter holds for the local gate family.
%   u=0 yields local gate i*I4 and an exactly known global phase.
%   n<=10 caps dense matrix storage; use factorized circuits for more.
%   References: J.-S. Caux, The Bethe Ansatz (trigonometric R-matrix);
%   E. Vernier et al., Phys. Rev. Lett. 133 (2024), 050606.
properties={'square','unitary','normal','complex','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'<=',20,'finite'});
if nargin<2 || isempty(u),u=0.6;end
if nargin<3 || isempty(gamma),gamma=pi/3;end
if nargin<4 || isempty(boundary)
    if n>=4 && mod(n,2)==0,boundary='periodic';else,boundary='open';end
end
validateattributes(u,{'numeric'},{'scalar','real','finite'});
validateattributes(gamma,{'numeric'},{'scalar','real','finite','>',0,'<',pi});
boundary=lower(char(boundary));
if ~any(strcmp(boundary,{'open','periodic'}))
    error('xxz_floquet:Boundary','Use open or periodic boundaries.');
end
if strcmp(boundary,'periodic') && (n<4 || mod(n,2)~=0)
    error('xxz_floquet:Boundary','Periodic brickwork requires even n>=4.');
end
if nargin<5,output='matrix';end
implicit=um_is_operator(output);
if ~implicit && n>10
    error('xxz_floquet:ExplicitOrder','Explicit matrix supports n<=10.');
end
gate=xxz_r_gate(u,gamma);
N=2^n;
gates=struct('gate',{},'sites',{});
for k=1:2:n-1
    gates(end+1)=struct('gate',gate,'sites',[k,k+1]); %#ok<AGROW>
end
for k=2:2:n-1
    gates(end+1)=struct('gate',gate,'sites',[k,k+1]); %#ok<AGROW>
end
if strcmp(boundary,'periodic')
    gates(end+1)=struct('gate',gate,'sites',[n,1]);
end
op=@(X,flag)um_circuit_apply(X,flag,n,gates,1);
if implicit,U=op;else,U=op(eye(N),'notransp');end
end
