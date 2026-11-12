function [U,properties] = cayley_jacobi(n,a,b,t)
%CAYLEY_JACOBI Cayley image of the symmetric Jacobi recurrence matrix.
%   Given the Jacobi probability weight (1-x)^A(1+x)^B on [-1,1],
%   A,B>-1, construct its N-by-N Jacobi recurrence J for x in [-2,2].
%       U=(J-1i*T*I)/(J+1i*T*I), T>0.
%   The eigenvalues are (x_j-1i*T)/(x_j+1i*T), x_j=eig(J).
%   This is a unitary image of an Hermitian matrix, not a finite CMV
%   or Gauss-Szego quadrature matrix.
%   Reference: Aurentz, Mach, Vandebril, Watkins, ETNA 46 (2017),
%   447-459 (Cayley eigenvalue methods).
properties={'square','unitary','normal','complex','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(a),a=0;end
if nargin<3 || isempty(b),b=0;end
if nargin<4 || isempty(t),t=1;end
validateattributes(t,{'numeric'},{'scalar','real','positive','finite'});
[d,s]=js_jacobi_recurrence(a,b,n);
J=diag(d)+diag(s,1)+diag(s,-1);
I=eye(n);
U=(J-1i*t*I)/(J+1i*t*I);
end
