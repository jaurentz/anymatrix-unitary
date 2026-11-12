function [U,properties] = grover_operator(n,marked)
%GROVER_OPERATOR Finite Grover search step, product of two reflections.
%   U=GROVER_OPERATOR(N,M) marks the first M basis vectors, 1<=M<N.
%   U=(2*s*s'-I)*(I-2*Pi), s=ones(N,1)/sqrt(N).
%   Exact spectrum: exp(+-2i*asin(sqrt(M/N))) once each, +1 of
%   multiplicity M-1, -1 of multiplicity N-M-1. Returns a dense matrix.
%   Reference: L. K. Grover, STOC (1996), pp. 212-219.
%   https://doi.org/10.1145/237814.237866
properties={'square','unitary','normal','real','orthogonal', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2 || isempty(marked),marked=1;end
validateattributes(marked,{'numeric'},{'scalar','integer','>=',1,'<',n});
s=ones(n,1)/sqrt(n);
reflection=2*(s*s')-eye(n);
oracle=diag([ -ones(marked,1); ones(n-marked,1)]);
U=reflection*oracle;
end
