function [U,properties] = finite_zak(n,p,output)
%FINITE_ZAK Finite, normalized Zak transform with rectangular indexing.
%   N=P*Q, with index j+P*k for j=0:P-1 and k=0:Q-1.
%   U applies a unitary length-Q DFT in the k coordinate at fixed j:
%   U=kron(F_Q,I_P). P defaults to 1; P=N gives identity.
%   Suitable for finite time-frequency/Gabor analysis.
%   Reference: H. Bolcskei and F. Hlawatsch, IEEE TSP 45 (1997).
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [51] https://doi.org/10.1109/78.564174
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(p),p=1;end
if nargin<3,output='matrix';end
if nargin==2 && (ischar(p)||isstring(p)),output=p;p=1;end
validateattributes(p,{'numeric'},{'scalar','integer','positive','finite'});
if mod(n,p)~=0
    error('finite_zak:Factor','P must divide N.');
end
q=n/p;
if um_is_operator(output)
    U=@(X,flag)um_zak_apply(X,flag,p,q);
    return
end
k=(0:q-1)';
F=exp(-2i*pi*(k*k')/q)/sqrt(q);
U=kron(F,eye(p));
end
