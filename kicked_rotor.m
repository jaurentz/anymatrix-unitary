function [U,properties] = kicked_rotor(n,K,hbar,theta_shift,output)
%KICKED_ROTOR Finite Fourier-grid Floquet quantum kicked rotor.
%   U = Dp*F'*Dx*F, with F(j,k)=exp(-2*pi*1i*j*k/n)/sqrt(n),
%   Dp(k,k)=exp(-1i*hbar*k^2/2), k=0,...,n-1,
%   Dx(j,j)=exp(-1i*K*cos(2*pi*(j+theta_shift)/n)/hbar).
%   K is classical kick strength; hbar is effective Planck's constant.
%   Default hbar=2*pi/n and theta_shift=0. This fixes a finite-grid
%   convention and is NOT an arbitrary truncation of an infinite operator.
%   Momentum labels are uncentered (0,...,n-1); different torus twists
%   and boundary conventions generally change the finite spectrum.
%   K=0 gives an exactly known diagonal spectrum.
%   Reference: F. M. Izrailev, Phys. Rep. 196 (1990), 299-392.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [35] https://doi.org/10.1016/0370-1573(90)90067-C
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'finite'});
if nargin<2 || isempty(K),K=5;end
if nargin<3 || isempty(hbar),hbar=2*pi/n;end
if nargin<4 || isempty(theta_shift),theta_shift=0;end
validateattributes(K,{'numeric'},{'scalar','real','finite'});
validateattributes(hbar,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(theta_shift,{'numeric'},{'scalar','real','finite'});
q=(0:n-1)';
Dp=exp(-1i*(hbar/2)*(q.^2));
Dx=exp(-1i*(K/hbar)*cos(2*pi*(q+theta_shift)/n));
if nargin<5,output='matrix';end
op=@(X,flag)um_fft_apply(X,flag,'kick',Dx,Dp);
if um_is_operator(output),U=op;else,U=op(eye(n),'notransp');end
end
