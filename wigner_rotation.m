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
%   For beta=k*pi with integer k, Euler's middle factor is EXACTLY
%   a monomial matrix in this m=-ell:ell order, for integer AND
%   half-integer spin:
%     beta=0 mod 4*pi: I; beta=2*pi mod 4*pi: (-1)^(2*ell)*I;
%     beta=pi mod 4*pi: S; beta=3*pi mod 4*pi: (-1)^(2*ell)*S;
%     S(row m,col m')=(-1)^(ell-m') delta_{m,-m'}.
%   Consequently special-angle full D is diagonal or anti-diagonal,
%   with O(ell) storage/action (up to FFT-free complex phases).
%   OUTPUT='operator' enables exact algebraic forward/adjoint/transpose
%   multi-RHS action for these beta ONLY, and rejects generic beta.
%   For ordinary 'matrix' at special beta, sparse output is returned
%   above order 512; at smaller order legacy dense output is retained.
%   The branch detects only EXACT represented beta=k*pi, not a
%   near-pi tolerance that would silently change the requested matrix.
%   Generic middle angles still use dense O((2*ELL+1)^3) expm setup;
%   no fast arbitrary-angle Wigner method is claimed.
%   Reference: Gimbutas and Greengard, JCP 228 (2009), 5621-5627.
%   [108] https://doi.org/10.1016/j.jcp.2009.05.014
%   [112] https://doi.org/10.1007/s00041-008-9013-5
%   Integer spin factors through SO(3); half-integer spin represents SU(2).
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(ell,{'numeric'}, ...
    {'scalar','real','nonnegative','finite'});
if 2*ell~=round(2*ell)
    error('wigner_rotation:Spin','Spin must be a nonnegative integer or half-integer.');
end
if nargin<2||isempty(alpha),alpha=.3;end
if nargin<3||isempty(beta),beta=.7;end
if nargin<4||isempty(gamma),gamma=.4;end
validateattributes(alpha,{'numeric'},{'scalar','real','finite'});
validateattributes(beta,{'numeric'},{'scalar','real','finite'});
validateattributes(gamma,{'numeric'},{'scalar','real','finite'});
if nargin<5||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('wigner_rotation:Mode','OUTPUT must be scalar text.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('wigner_rotation',ell,alpha,beta,gamma);
    return
end
if ~any(strcmp(output,{'matrix','operator'}))
    error('wigner_rotation:Mode', ...
        'OUTPUT must be matrix, operator, or spectral.');
end
n=round(2*ell)+1;
m=(-ell:ell)';
if ell==0
    if strcmp(output,'operator')
        U=@(X,flag)um_wigner_axis_apply(X,flag,1,1);
    else
        U=1;
    end
    return
end
q=round(beta/pi);
axisExact=abs(q)<=2^20 && beta==q*pi;
if axisExact
    eta=(-1)^round(2*ell); % 2*pi SU(2) phase for half-integer spins
    r=mod(q,4);
    if mod(r,2)==0
        perm=(1:n)';
        spinFactor=1;
        if r==2,spinFactor=eta;end
        phase=spinFactor*exp(-1i*(alpha+gamma)*m);
    else
        perm=(n:-1:1)';
        mCol=m(perm);
        spinFactor=1;
        if r==3,spinFactor=eta;end
        phase=spinFactor*((-1).^round(ell-mCol)) ...
            .*exp(-1i*(alpha*m+gamma*mCol));
    end
    if strcmp(output,'operator')
        U=@(X,flag)um_wigner_axis_apply(X,flag,phase,perm);
    else
        U=sparse((1:n)',perm,phase,n,n);
        if n<=512,U=full(U);end
    end
    return
end
if strcmp(output,'operator')
    error('wigner_rotation:Operator', ...
        'Fast operator is available only for exactly represented beta=k*pi.');
end
c=sqrt((ell-m(1:end-1)).*(ell+m(1:end-1)+1));
K=(diag(c,1)-diag(c,-1))/2;
U=diag(exp(-1i*alpha*m))*expm(beta*K)*diag(exp(-1i*gamma*m));
end
