function [U,properties] = selfdual_kicked_ising(n,fields)
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
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',3,'<=',11,'finite'});
if nargin<2 || isempty(fields),fields=0.3*sin(1:n);end
validateattributes(fields,{'numeric'},{'real','finite','vector'});
if isscalar(fields),fields=repmat(fields,1,n);end
if numel(fields)~=n
    error('selfdual_kicked_ising:Fields','Fields must be scalar or length n.');
end
fields=double(fields(:));
N=2^n;
states=uint64((0:N-1)');
z=zeros(N,n);
for j=1:n
    z(:,j)=1-2*double(bitget(states,j));
end
energy=sum(z.*z(:,[2:n,1]),2)+z*fields;
phase=exp(-1i*(pi/4)*energy);
R=[1,-1i;-1i,1]/sqrt(2);
Ux=1;
for j=1:n
    Ux=kron(Ux,R);
end
U=bsxfun(@times,Ux,phase.');
end
