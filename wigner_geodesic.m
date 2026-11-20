function [U,properties]=wigner_geodesic(n,alpha,beta,gamma,output)
%WIGNER_GEODESIC Wigner D representation of rotations on SU(2)/SO(3).
%   U=WIGNER_GEODESIC(N,ALPHA,BETA,GAMMA) constructs the spin-j
%   representation of the Z-Y-Z Euler rotation, j=(N-1)/2. The
%   representation of the compact rotation group is unitary. Every
%   element can be reached by a geodesic from identity in its usual
%   bi-invariant metric; this generic Euler coordinate path itself is
%   not a geodesic. Half-integer j (even N) represents SU(2), not SO(3).
%   U=exp(-i ALPHA Jz)*exp(-i BETA Jy)*exp(-i GAMMA Jz).
%   The middle factor is an explicit real orthogonal matrix obtained
%   by exponentiating the sparse angular-momentum tridiagonal generator.
%   Exact eigenvalues exp(-i*m*theta), m=-j,...,j, where
%   cos(theta/2)=cos(beta/2)*cos((alpha+gamma)/2),
%   theta in [0,2pi] chosen to retain the SU(2) lift.
%   This is dense by default; no fast matrix-free mode is asserted.
%   [99] Kostelec and Rockmore, J. Fourier Anal. Appl. 14 (2008).
%   https://doi.org/10.1007/s00041-008-9013-5
properties={'square','unitary','normal','complex','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(alpha),alpha=.31;end
if nargin<3||isempty(beta),beta=.77;end
if nargin<4||isempty(gamma),gamma=-.46;end
for a=[alpha,beta,gamma]
    validateattributes(a,{'numeric'},{'scalar','real','finite'});
end
if nargin<5||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('wigner_geodesic:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
j=(n-1)/2;m=(-j:j)';
if strcmp(output,'spectral')
    c=max(-1,min(1,cos(beta/2)*cos((alpha+gamma)/2)));
    theta=2*acos(c);
    z=exp(-1i*m*theta);
    U=unitary_spectral_reference('wigner_geodesic',n);
    U.finite.kind='exact';U.finite.eigenvalues=@()z;
    U.finite.moments=@(k)reshape(arrayfun(@(ell) ...
        mean(z.^ell),k(:)),size(k));
    U.reference='Kostelec-Rockmore, JFAA 14 (2008), doi:10.1007/s00041-008-9013-5';
    U.notes='Exact spin-j character/eigenphase law; even N uses SU(2) double cover.';
elseif any(strcmp(output,{'matrix','explicit'}))
    if n==1
        U=1;return
    end
    coeff=sqrt((j-m(1:end-1)).*(j+m(1:end-1)+1));
    P=sparse(2:n,1:n-1,coeff,n,n);
    B=expm(full((beta/2)*(P'-P)));
    U=exp(-1i*alpha*m).*B.*(exp(-1i*gamma*m)).';
else
    error('wigner_geodesic:Mode','OUTPUT must be matrix or spectral.');
end
end
