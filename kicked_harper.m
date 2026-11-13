function [U,properties] = kicked_harper(n,K,L,output)
%KICKED_HARPER Finite torus version of the kicked Harper Floquet operator.
%   U=D_p(L)*F'*D_x(K)*F, with D_x and D_p diagonal phase matrices
%   exp(-1i*K*cos(2*pi*k/N)) and exp(-1i*L*cos(2*pi*k/N)).
%   F is the normalized negative-exponent DFT (N-by-N).
%   This exact finite unitary definition fixes phase and ordering
%   conventions; its spectral behavior varies with K,L and dimension.
%   Reference: Artuso, Casati, Shepelyansky, Phys. Rev. Lett. 68
%   (1992), 3826-3829. https://doi.org/10.1103/PhysRevLett.68.3826
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=4 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('kicked_harper',n,K,L);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(K),K=1.2;end
if nargin<3 || isempty(L),L=0.7;end
validateattributes(K,{'numeric'},{'scalar','real','finite'});
validateattributes(L,{'numeric'},{'scalar','real','finite'});
q=(0:n-1)';
Dx=exp(-1i*K*cos(2*pi*q/n));
Dp=exp(-1i*L*cos(2*pi*q/n));
if nargin<4,output='matrix';end
op=@(X,flag)um_fft_apply(X,flag,'kick',Dx,Dp);
if um_is_operator(output),U=op;else,U=op(eye(n),'notransp');end
end
