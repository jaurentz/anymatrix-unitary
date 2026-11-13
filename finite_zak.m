function [U,properties] = finite_zak(n,p)
%FINITE_ZAK Finite, normalized Zak transform with rectangular indexing.
%   N=P*Q, with index j+P*k for j=0:P-1 and k=0:Q-1.
%   U applies a unitary length-Q DFT in the k coordinate at fixed j:
%   U=kron(F_Q,I_P). P defaults to 1; P=N gives identity.
%   Suitable for finite time-frequency/Gabor analysis.
%   Reference: H. Bolcskei and F. Hlawatsch, IEEE TSP 45 (1997).
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(p),p=1;end
validateattributes(p,{'numeric'},{'scalar','integer','positive','finite'});
if mod(n,p)~=0
    error('finite_zak:Factor','P must divide N.');
end
q=n/p;
k=(0:q-1)';
F=exp(-2i*pi*(k*k')/q)/sqrt(q);
U=kron(F,eye(p));
end
