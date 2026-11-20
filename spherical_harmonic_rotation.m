function [U,properties]=spherical_harmonic_rotation(L,alpha,beta,gamma,output)
%SPHERICAL_HARMONIC_ROTATION Orthogonal direct sum of Wigner blocks.
%   U=SPHERICAL_HARMONIC_ROTATION(L,ALPHA,BETA,GAMMA) represents a
%   single SO(3) rotation on all complex spherical harmonics of
%   degrees ell=0,...,L-1, for the z-y-z Euler-angle convention.
%   Its order is L^2. The mathematical matrix is unitary; in the
%   complex harmonic basis it is generally not real.
%   With physical rotation angle delta, exp(-i*m*delta) occurs
%   with multiplicity L-|m|, m=-(L-1),...,L-1. Hence
%       tr(U^k)=|sum(j=0:L-1,exp(-i*k*j*delta))|^2.
%   In the special regime delta=c/L, 0<c<pi, the empirical
%   eigenphase density tends to (1-|theta|/c)/c on |theta|<=c.
%   Sparse block-diagonal explicit representation; no fast
%   spherical transform (FFT/butterfly) is claimed.
%   [97] https://doi.org/10.1016/j.jcp.2009.05.014
properties={'square','unitary','normal','sparse', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(alpha),alpha=.3;end
if nargin<3||isempty(beta),beta=.7;end
if nargin<4||isempty(gamma),gamma=.4;end
validateattributes(alpha,{'numeric'},{'scalar','real','finite'});
validateattributes(beta,{'numeric'},{'scalar','real','finite'});
validateattributes(gamma,{'numeric'},{'scalar','real','finite'});
if nargin>=5
    if ~(ischar(output)||isstring(output)) || ~strcmpi(char(output),'spectral')
        error('spherical_harmonic_rotation:Mode','Only terminal spectral mode is supported.');
    end
    U=unitary_spectral_reference('spherical_harmonic_rotation',L,alpha,beta,gamma);
    return
end
N=L^2;
nStored=sum((2*(0:L-1)+1).^2);
U=spalloc(N,N,nStored);
offset=0;
for ell=0:L-1
    idx=offset+(1:2*ell+1);
    U(idx,idx)=sparse(wigner_rotation(ell,alpha,beta,gamma));
    offset=offset+2*ell+1;
end
end
