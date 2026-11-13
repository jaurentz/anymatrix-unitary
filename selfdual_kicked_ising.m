function [U,properties] = selfdual_kicked_ising(n,fields,output)
%SELFDUAL_KICKED_ISING Self-dual periodically kicked Ising spin chain.
%   U = exp(-1i*pi/4*sum(X_j)) *
%       exp(-1i*pi/4*(sum(Z_j Z_(j+1))+sum(fields_j*Z_j))).
%   Periodic boundaries. n spin sites, dimension 2^n (n>=3).
%   fields is a real scalar or a real n-vector of longitudinal fields.
%   Default fields_j=0.3*sin(j), deterministic and inhomogeneous.
%   Basis: bit j-1 corresponds to spin j; Z|0>=|0>.
%   The model is self-dual in the coupling/kick strengths pi/4;
%   the finite-size random-matrix form-factor claims require additional
%   symmetry, disorder, and limiting assumptions.
%   Reference: B. Bertini, P. Kos, T. Prosen, Phys. Rev. Lett. 121
%   (2018), 264101. https://doi.org/10.1103/PhysRevLett.121.264101
properties={'square','unitary','normal','complex','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('selfdual_kicked_ising',n,fields);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',3,'<=',20,'finite'});
if nargin<2 || isempty(fields),fields=0.3*sin(1:n);end
validateattributes(fields,{'numeric'},{'real','finite','vector'});
if isscalar(fields),fields=repmat(fields,1,n);end
if numel(fields)~=n
    error('selfdual_kicked_ising:Fields','Fields must be scalar or length n.');
end
fields=double(fields(:));
if nargin<3,output='matrix';end
implicit=um_is_operator(output);
if ~implicit && n>11
    error('selfdual_kicked_ising:ExplicitOrder','Explicit matrix supports n<=11.');
end
N=2^n;
states=uint64((0:N-1)');
z=zeros(N,n);
for j=1:n,z(:,j)=1-2*double(bitget(states,j));end
energy=sum(z.*z(:,[2:n,1]),2)+z*fields;
phase=exp(-1i*(pi/4)*energy);
R=[1,-1i;-1i,1]/sqrt(2);
op=@(X,flag)um_spin_apply(X,flag,n,phase,R,'diag_right');
if implicit,U=op;else,U=op(eye(N),'notransp');end
end
