function [U,properties] = entangled_reflection(d,phi,output)
%ENTANGLED_REFLECTION Bipartite maximally entangled phase reflection.
%   v=vec(I_d)/sqrt(d), U=I+(exp(1i*phi)-1)*v*v' has order d^2.
%   Default PHI=pi gives the reflection from Kukulski et al. (2024).
%   Spectrum: exp(1i*phi) once, 1 with multiplicity d^2-1.
%   At PHI=pi and d>=2, the product numerical range is [1-2/d,1].
%   For d>=2, ordinary numerical range contains zero at PHI=pi.
%   OUTPUT='operator' returns a function supporting ('size'),
%   ('notransp'), and ('transp' or 'ctransp') with multiple RHS.
%   Ref: R. Kukulski et al., Sci. Rep. 14 (2024), 26588, Eq. (15).
%   https://doi.org/10.1038/s41598-024-75148-z
properties={'square','unitary','normal','sparse','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(phi),phi=pi;end
validateattributes(phi,{'numeric'},{'scalar','real','finite'});
if nargin<3||isempty(output),output='matrix';end
n=d*d;
if ~isfinite(n)||n>flintmax
    error('entangled_reflection:Order','d^2 exceeds the exact index range.');
end
ii=1+(0:d-1)'*(d+1);
v=sparse(ii,ones(d,1),ones(d,1)/sqrt(d),n,1);
if phi==pi,lambda=-1;else,lambda=exp(1i*phi);end
if um_is_operator(output)
    U=@(X,flag)um_entangled_reflection_apply(X,flag,v,lambda,n);
else
    U=speye(n)+(lambda-1)*(v*v');
end
end

