function [U,properties]=real_spherical_harmonic_rotation(L,alpha,beta,gamma,output)
%REAL_SPHERICAL_HARMONIC_ROTATION Real orthogonal SO(3) harmonics.
%   U=REAL_SPHERICAL_HARMONIC_ROTATION(L,A,B,G) represents z-y-z
%   SO(3) Euler rotation on ALL *real* spherical harmonics of degree
%   ell=0,...,L-1. U has order L^2 and is real orthogonal, not the
%   entrywise real part of the complex Wigner representation.
%   For each ell, Dreal=C_ell'*Dcomplex*C_ell with the real
%   Condon--Shortley sine/cosine basis in UM_WIGNER_REAL_BASIS.
%   Degree block order m=-ell:ell has sine for m<0, constant
%   for m=0, cosine for m>0. This basis convention fixes signs.
%
%   OUTPUT='matrix' (default): sparse block-diagonal real orthogonal
%       matrix, explicit O(L^3) stored entries and O(L^4) dense
%       block-exponential setup (not FFT-fast construction).
%       Dense assembly limited to L<=72 to bound cubic storage.
%   OUTPUT='operator': sparse real-to-complex Condon--Shortley
%       basis conversions around matrix-free skew-tridiagonal Wigner
%       exponential actions, with no precomputed dense degree blocks.
%       For fixed beta this is O(L^3) work per RHS, O(L^2) stored
%       data and O(L^2) setup, compared with O(L^4) dense expm
%       setup in 'matrix'. The action is approximate in floating
%       point and NOT a subquadratic single-degree Wigner FFT.
%       Real inputs yield real outputs (up to discarded roundoff),
%       complex multiple RHS remain complex-linear.
%   OUTPUT='spectral': exact finite phase multiset/trace moments
%       inherited from the equivalent complex SO(3) representation.
%       At physical angle delta, e^(-im*delta) has multiplicity L-|m|.
%       A real orthogonal representation has conjugate phase pairs.
%   OUTPUT='basis': unitary Condon--Shortley basis changes C_ell
%       for each degree, useful for verifying convention and blocks.
%
%   [108] Gimbutas and Greengard, JCP 228 (2009), 5621-5627.
%   https://doi.org/10.1016/j.jcp.2009.05.014
%   [112] Real harmonic representation background:
%   https://doi.org/10.1007/s00041-008-9013-5
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','sparse','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',128});
if nargin<2||isempty(alpha),alpha=.3;end
if nargin<3||isempty(beta),beta=.7;end
if nargin<4||isempty(gamma),gamma=.4;end
validateattributes(alpha,{'numeric'},{'scalar','real','finite'});
validateattributes(beta,{'numeric'},{'scalar','real','finite'});
validateattributes(gamma,{'numeric'},{'scalar','real','finite'});
if nargin<5||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('real_spherical_harmonic_rotation:Mode', ...
        'OUTPUT must be a scalar string or character row.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','operator','spectral','basis'}))
    error('real_spherical_harmonic_rotation:Mode','Unknown output mode.');
end
if strcmp(output,'spectral')
    U=unitary_spectral_reference('spherical_harmonic_rotation', ...
        L,alpha,beta,gamma);
    U.family='real_spherical_harmonic_rotation';
    U.reference='Real Condon--Shortley basis of spherical Wigner D (JCP 2009)';
    U.notes=[U.notes ' Real orthogonal SO(3) conjugate-paired ', ...
        'representation; unchanged trace eigenphases.'];
    return
end
if strcmp(output,'matrix')&&L>72
    error('real_spherical_harmonic_rotation:ExplicitSize', ...
        'Use operator or spectral above L=72.');
end
if strcmp(output,'basis')
    U=cell(L,1);
    for ell=0:L-1,U{ell+1}=um_wigner_real_basis(ell);end
    return
end
if strcmp(output,'operator')
    blocks=cell(L,1);
    for ell=0:L-1
        C=um_wigner_real_basis(ell);
        D=wigner_rotation(ell,alpha,beta,gamma,'operator');
        blocks{ell+1}=@(X,flag)um_wigner_real_block_apply(X,flag,C,D);
    end
    U=@(X,flag)um_wigner_direct_sum_apply(X,flag,blocks,L^2);
    return
end
blocks=cell(L,1);
for ell=0:L-1
    C=um_wigner_real_basis(ell);
    Dc=wigner_rotation(ell,alpha,beta,gamma);
    Dr=C'*Dc*C;
    if norm(imag(Dr),'fro')>2e-10*(2*ell+1)
        error('real_spherical_harmonic_rotation:RealBasis', ...
            'Condon--Shortley real-basis conjugate symmetry failed.');
    end
    blocks{ell+1}=real(Dr);
end
N=L^2;
nz=sum((2*(0:L-1)+1).^2);
U=spalloc(N,N,nz);
offset=0;
for ell=0:L-1
    idx=offset+(1:2*ell+1);
    U(idx,idx)=sparse(blocks{ell+1});
    offset=offset+2*ell+1;
end
end
