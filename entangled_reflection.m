function [U,properties] = entangled_reflection(d,phi,output)
%ENTANGLED_REFLECTION Maximally entangled bipartite phase reflection.
%   v=vec(I_d)/sqrt(d); U=I+(exp(1i*phi)-1)*v*v' has order d^2.
%   The default phi=pi yields I-2vv' with a product-range gap.
%   Spectrum exp(i*phi) once, 1 with multiplicity d^2-1.
%   OUTPUT='operator' supports forward/adjoint multi-RHS actions;
%   OUTPUT='spectral' returns exact finite eigenphase data.
%   Ref: R. Kukulski et al., Scientific Reports 14 (2024).
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
if (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('entangled_reflection',d,phi);
    return
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
