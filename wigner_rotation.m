function [U,properties]=wigner_rotation(ell,alpha,beta,gamma,output)
%WIGNER_ROTATION Unitary degree-ell spherical-harmonic SO(3) rotation.
%   U=WIGNER_ROTATION(ELL,ALPHA,BETA,GAMMA) returns D^ell(R) for
%   integer ELL>=0 and z-y-z Euler angles (radians). The ordered
%   complex spherical-harmonic basis has m=-ELL,...,ELL.
%   U=exp(-i*alpha*Jz)*exp(-i*beta*Jy)*exp(-i*gamma*Jz).
%   The middle factor is the matrix exponential of the explicit
%   real skew-symmetric angular-momentum generator (no Wigner
%   factorial sum). The result is mathematically unitary.
%   Every SO(3) rotation has an angle delta, giving exact eigenvalues
%   exp(-i*m*delta) for m=-ELL:ELL; for this z-y-z convention,
%   cos(delta/2)=cos(beta/2)*cos((alpha+gamma)/2).
%   Defaults ALPHA=.3, BETA=.7, GAMMA=.4.
%   This is a dense O((2*ELL+1)^3) setup, not a fast rotation kernel.
%   Reference: Gimbutas and Greengard, JCP 228 (2009), 5621-5627.
%   [108] https://doi.org/10.1016/j.jcp.2009.05.014
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(ell,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
if nargin<2||isempty(alpha),alpha=.3;end
if nargin<3||isempty(beta),beta=.7;end
if nargin<4||isempty(gamma),gamma=.4;end
validateattributes(alpha,{'numeric'},{'scalar','real','finite'});
validateattributes(beta,{'numeric'},{'scalar','real','finite'});
validateattributes(gamma,{'numeric'},{'scalar','real','finite'});
if nargin>=5
    if ~(ischar(output)||isstring(output)) || ~strcmpi(char(output),'spectral')
        error('wigner_rotation:Mode','Only terminal spectral mode is supported.');
    end
    U=unitary_spectral_reference('wigner_rotation',ell,alpha,beta,gamma);
    return
end
if ell==0,U=1;return;end
m=(-ell:ell)';
c=sqrt((ell-m(1:end-1)).*(ell+m(1:end-1)+1));
K=(diag(c,1)-diag(c,-1))/2;
U=diag(exp(-1i*alpha*m))*expm(beta*K)*diag(exp(-1i*gamma*m));
end
