function [U,properties] = floquet_ising(n,J,h,boundary,output)
%FLOQUET_ISING Integrable periodically kicked transverse-field Ising chain.
%   U = exp(-1i*J*sum(Z_j Z_(j+1)))*exp(-1i*h*sum(X_j)).
%   n denotes spin sites; order(U)=2^n. j=n couples to j=1 for
%   'periodic' boundaries (default), otherwise use 'open'.
%   Computational basis uses bit j-1 for spin j; Z|0>=|0>.
%   The left factor is a diagonal Ising interaction layer and the
%   right factor is a Kronecker product of single-spin X rotations.
%   Dense assembly scales as 4^n; cap n at 11 for practical memory.
%   References: E. Lieb, T. Schultz, D. Mattis, Ann. Phys. 16 (1961);
%   P. Pfeuty, Ann. Phys. 57 (1970), 79-90.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [30] https://doi.org/10.1016/0003-4916(61)90115-4
%   [31] https://doi.org/10.1016/0003-4916(70)90270-8
properties={'square','unitary','normal','complex','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=5 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('floquet_ising',n,J,h,boundary);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'<=',20,'finite'});
if nargin<2 || isempty(J),J=0.65;end
if nargin<3 || isempty(h),h=0.43;end
if nargin<4 || isempty(boundary),boundary='periodic';end
validateattributes(J,{'numeric'},{'scalar','real','finite'});
validateattributes(h,{'numeric'},{'scalar','real','finite'});
boundary=lower(char(boundary));
if ~any(strcmp(boundary,{'open','periodic'}))
    error('floquet_ising:Boundary','Boundary must be open or periodic.');
end
if nargin<5,output='matrix';end
implicit=um_is_operator(output);
if ~implicit && n>11
    error('floquet_ising:ExplicitOrder','Explicit matrix supports n<=11.');
end
N=2^n;
states=uint64((0:N-1)');
z=zeros(N,n);
for j=1:n,z(:,j)=1-2*double(bitget(states,j));end
energy=sum(z(:,1:n-1).*z(:,2:n),2);
if strcmp(boundary,'periodic'),energy=energy+z(:,n).*z(:,1);end
phase=exp(-1i*J*energy);
R=[cos(h),-1i*sin(h);-1i*sin(h),cos(h)];
op=@(X,flag)um_spin_apply(X,flag,n,phase,R,'diag_left');
if implicit,U=op;else,U=op(eye(N),'notransp');end
end
