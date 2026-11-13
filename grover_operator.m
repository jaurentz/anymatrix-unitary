function [U,properties] = grover_operator(n,marked,output)
%GROVER_OPERATOR Finite Grover search step, product of two reflections.
%   U=GROVER_OPERATOR(N,M) marks the first M basis vectors, 1<=M<N.
%   U=(2*s*s'-I)*(I-2*Pi), s=ones(N,1)/sqrt(N).
%   Exact spectrum: exp(+-2i*asin(sqrt(M/N))) once each, +1 of
%   multiplicity M-1, -1 of multiplicity N-M-1. Returns a dense matrix.
%   Reference: L. K. Grover, STOC (1996), pp. 212-219.
%   https://doi.org/10.1145/237814.237866
properties={'square','unitary','normal','real','orthogonal', ...
    'scalable','eigenvalues','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('grover_operator',n,marked);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2 || isempty(marked),marked=1;end
validateattributes(marked,{'numeric'},{'scalar','integer','>=',1,'<',n});
if nargin<3,output='matrix';end
op=@(X,flag)um_reflection_apply(X,flag,n,marked);
if um_is_operator(output),U=op;else,U=op(eye(n),'notransp');end
end
