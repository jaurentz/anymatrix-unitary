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
%   'operator' uses NO dense degree blocks: each Wigner D is applied
%   by an O(ell) storage degree-32 scaled Taylor action of its real
%   skew tridiagonal generator. Fixed reduced beta costs O(ell^2)
%   per degree and O(L^3) for all L blocks, O(L^2) total stored
%   angular-momentum data, WITHOUT O(L^4) dense expm setup.
%   It is approximate in floating point, NOT a Gimbutas-Greengard
%   O(L^2 log L) Wigner transform. Transpose and ADJOINT flags differ.
%   For a real-orthogonal counterpart use the separate generator
%   REAL_SPHERICAL_HARMONIC_ROTATION with explicit Condon-Shortley
%   complex-to-real basis changes (not entrywise real(D)).
%   [108] https://doi.org/10.1016/j.jcp.2009.05.014
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
if nargin>=5 && ~isempty(output)
    if ~(ischar(output)||(isstring(output)&&isscalar(output)))
        error('spherical_harmonic_rotation:Mode','Invalid output.');
    end
    switch lower(char(output))
        case 'spectral'
            U=unitary_spectral_reference('spherical_harmonic_rotation', ...
                L,alpha,beta,gamma);
            return
        case 'operator'
            blocks=cell(L,1);
            for ell=0:L-1
                blocks{ell+1}=wigner_rotation(ell,alpha,beta,gamma,'operator');
            end
            U=@(X,flag)um_wigner_direct_sum_apply(X,flag,blocks,L^2);
            return
        case 'matrix'
        otherwise
            error('spherical_harmonic_rotation:Mode', ...
                'Use matrix, operator or spectral.');
    end
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
