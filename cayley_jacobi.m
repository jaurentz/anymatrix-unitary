function [U,properties] = cayley_jacobi(n,a,b,t,output)
%CAYLEY_JACOBI Cayley image of the symmetric Jacobi recurrence matrix.
%   Given the Jacobi probability weight (1-x)^A(1+x)^B on [-1,1],
%   A,B>-1, construct its N-by-N Jacobi recurrence J for x in [-2,2].
%       U=(J-1i*T*I)/(J+1i*T*I), T>0.
%   The eigenvalues are (x_j-1i*T)/(x_j+1i*T), x_j=eig(J).
%   This is a unitary image of an Hermitian matrix, not a finite CMV
%   or Gauss-Szego quadrature matrix.
%   Reference: Aurentz, Mach, Vandebril, Watkins, ETNA 46 (2017),
%   447-459 (Cayley eigenvalue methods).
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [23] https://www.maths.tcd.ie/EMIS/journals/ETNA/volumes/2011-2020/vol46/abstract_vol46_pp447-459.html?pages=447-459&vol=46
properties={'square','unitary','normal','complex','scalable', ...
    'eigenvalues','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=5 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('cayley_jacobi',n,a,b,t);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(a),a=0;end
if nargin<3 || isempty(b),b=0;end
if nargin<4 || isempty(t),t=1;end
validateattributes(t,{'numeric'},{'scalar','real','positive','finite'});
[d,s]=js_jacobi_recurrence(a,b,n);
if nargin<5,output='matrix';end
J=spdiags(d(:),0,n,n);
if n>1
    J=J+sparse(1:n-1,2:n,s(:),n,n)+sparse(2:n,1:n-1,s(:),n,n);
end
% Explicit index lists avoid ambiguity in spdiags off-diagonal alignment.
if um_is_operator(output)
    plusFactor=decomposition(J+1i*t*speye(n),'lu');
    minusFactor=decomposition(J-1i*t*speye(n),'lu');
    U=@(X,flag)um_cayley_apply(X,flag,plusFactor,minusFactor,t,n);
else
    I=eye(n);
    U=(full(J)-1i*t*I)/(full(J)+1i*t*I);
end
end
