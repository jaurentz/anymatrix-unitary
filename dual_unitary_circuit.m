function [U,properties] = dual_unitary_circuit(n,depth,chi,angles,boundary)
%DUAL_UNITARY_CIRCUIT Brickwork Floquet circuit of dual-unitary qubit gates.
%   [U,properties]=DUAL_UNITARY_CIRCUIT(n,depth,chi,angles,boundary).
%   n is spin count (2<=n<=10), U has order 2^n.
%   Each gate is (Ry(a) kron Rz(b))*SWAP*exp(-i*chi*ZZ)*
%   (Rz(c) kron Ry(d)), with ordinary Kronecker products.
%   Default boundary is periodic for even n>=4, open otherwise.
%   This gate and its 4-by-4 index realignment are both unitary.
%   Odd bonds (1,2),(3,4),... act first, then even bonds
%   (2,3),(4,5),... and (n,1) for periodic even n>=4.
%   The same two-layer Floquet unitary is repeated depth times.
%   Open boundaries support n>=2; periodic boundaries require even n>=4.
%   The gate family is dual-unitary, but no universal CUE spectral law
%   is asserted for a finite, spatially homogeneous circuit.
%   References: B. Bertini, P. Kos, T. Prosen, Phys. Rev. Lett.
%   123 (2019), 210601; Commun. Math. Phys. 387 (2021), 597-620.
properties={'square','unitary','normal','complex','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'<=',10,'finite'});
if nargin<2 || isempty(depth),depth=1;end
if nargin<3 || isempty(chi),chi=0.34;end
if nargin<4 || isempty(angles),angles=[0.37,-0.61,0.29,0.83];end
if nargin<5 || isempty(boundary)
    if n>=4 && mod(n,2)==0,boundary='periodic';else,boundary='open';end
end
validateattributes(depth,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(chi,{'numeric'},{'scalar','real','finite'});
validateattributes(angles,{'numeric'},{'vector','numel',4,'real','finite'});
boundary=lower(char(boundary));
if ~any(strcmp(boundary,{'open','periodic'}))
    error('dual_unitary_circuit:Boundary','Use open or periodic boundaries.');
end
if strcmp(boundary,'periodic') && (n<4 || mod(n,2)~=0)
    error('dual_unitary_circuit:Boundary', ...
          'Periodic brickwork requires even n>=4.');
end
N=2^n;
if depth==0,U=eye(N);return;end
a=angles(:);
Ry1=[cos(a(1)/2),-sin(a(1)/2);sin(a(1)/2),cos(a(1)/2)];
Ry2=[cos(a(4)/2),-sin(a(4)/2);sin(a(4)/2),cos(a(4)/2)];
Rz1=diag(exp(-1i*a(2)*[1;-1]/2));
Rz2=diag(exp(-1i*a(3)*[1;-1]/2));
Swap=[1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1];
D=diag(exp(-1i*chi*[1;-1;-1;1]));
gate=kron(Ry1,Rz1)*Swap*D*kron(Rz2,Ry2);
odd=speye(N);
for k=1:2:n-1
    odd=qm_two_site_gate(n,gate,k,k+1)*odd;
end
even=speye(N);
for k=2:2:n-1
    even=qm_two_site_gate(n,gate,k,k+1)*even;
end
if strcmp(boundary,'periodic')
    even=qm_two_site_gate(n,gate,n,1)*even;
end
U=full(even*odd)^depth;
end
