function [U,properties]=discrete_sine(n,kind,output)
%DISCRETE_SINE Orthonormal discrete sine transform, types I--IV.
%   DISCRETE_SINE(N,KIND) returns a real orthogonal N-by-N transform.
%   KIND defaults to 2. DST-III is the inverse/transpose of DST-II.
%   DISCRETE_SINE(N,KIND,'operator') applies the FFT-based actions.
%   Reference: M. Puschel and J. M. F. Moura, IEEE TSP 56 (2008).
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [54] https://doi.org/10.1016/b978-0-12-373624-6.x5000-0
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(kind),kind=2;end
if nargin<3,output='matrix';end
if nargin==2 && (ischar(kind)||isstring(kind)),output=kind;kind=2;end
validateattributes(kind,{'numeric'},{'scalar','integer','>=',1,'<=',4});
if um_is_operator(output)
    U=@(X,flag)um_dst_apply(X,flag,n,kind);
    return
end
[K,J]=ndgrid(0:n-1,0:n-1);
switch kind
    case 1
        U=sqrt(2/(n+1))*sin(pi*(K+1).*(J+1)/(n+1));
    case {2,3}
        c=ones(n,1);c(end)=1/sqrt(2);
        S=sqrt(2/n)*bsxfun(@times,c, ...
            sin(pi*(K+1).*(J+0.5)/n));
        if kind==2,U=S;else,U=S';end
    case 4
        U=sqrt(2/n)*sin(pi*(K+0.5).*(J+0.5)/n);
end
end
