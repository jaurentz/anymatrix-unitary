function [U,properties] = cayley_hermitian_pencil(n,a,b,t,mode)
%CAYLEY_HERMITIAN_PENCIL Strict unitary Cayley pencil from a Jacobi matrix.
%   For real symmetric Jacobi H on [-2,2], A=H-1i*T*I and
%   B=H+1i*T*I, T>0. Then B\A is unitary and AA'=BB'.
%   The output equals CAYLEY_JACOBI(N,A,B,T) but the actual pencil
%   coefficients are accessible as outputs three and four.
%   Default Jacobi exponents a=b=0, T=1.
%   Aurentz et al., ETNA 46 (2017), 447-459.
%   [23] https://www.maths.tcd.ie/EMIS/journals/ETNA/volumes/2011-2020/vol46/abstract_vol46_pp447-459.html?pages=447-459&vol=46
%   Terminal 'pencil' mode returns struct P with P.A, P.B, P.U.
%   Standard [U,PROPERTIES] signature preserves Anymatrix discovery.
properties={'square','unitary','normal','complex','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,U=[];A=[];B=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(a),a=0;end
if nargin<3 || isempty(b),b=0;end
if nargin<4 || isempty(t),t=1;end
validateattributes(a,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(b,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(t,{'numeric'},{'scalar','real','finite','positive'});
if nargin<5 || isempty(mode),mode='matrix';end
if strcmpi(char(mode),'spectral')
    U=ep_spectral_reference('cayley_hermitian_pencil',n,a,b,t);
    A=[];B=[];return
end
if ~any(strcmpi(char(mode),{'matrix','pencil'}))
    error('cayley_hermitian_pencil:Mode','Supported modes: matrix, pencil, spectral.');
end
[d,s]=js_jacobi_recurrence(a,b,n);
H=spdiags(d(:),0,n,n);
if n>1
    H=H+sparse(1:n-1,2:n,s(:),n,n)+sparse(2:n,1:n-1,s(:),n,n);
end
A=H-1i*t*speye(n);
B=H+1i*t*speye(n);
U=B\A;
if strcmpi(char(mode),'pencil')
    U=struct('A',A,'B',B,'U',U);
end
end
