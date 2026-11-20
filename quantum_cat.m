function [U,properties]=quantum_cat(n,s,output)
%QUANTUM_CAT Hannay-Berry quantization of a hyperbolic torus cat map.
%   U=QUANTUM_CAT(N,S) quantizes the classical symplectic map
%           T_S = [2*S, 1; 4*S-1, 2],   S=1,2,... .
%   T_S has determinant 1, trace 2*S+2>2, and the checkerboard parity
%   needed for periodic-position-and-momentum quantization for all N.
%   The default S=1 is T=[2 1;3 2], not the familiar [2 1;1 1],
%   which is incompatible with periodic boundary conditions for all N
%   in Hannay and Berry's convention.
%
%   For j,k=0:N-1, the kernel is
%      U(j+1,k+1) = exp(i*pi/4)/sqrt(N) ...
%                    * exp(2*pi*i*(j^2-j*k+S*k^2)/N).
%   This is a diagonal chirp, normalized inverse-sign DFT, and another
%   diagonal chirp. It is unitary for every positive N; all entries have
%   magnitude 1/sqrt(N). OUTPUT='operator' gives O(N log N) FFT actions.
%   The overall exp(i*pi/4) phase is part of this specific convention.
%   Generic exact finite eigenangle multiplicities are not asserted.
%
%   J. H. Hannay and M. V. Berry, Quantization of linear maps on a
%   torus: Fresnel diffraction by a periodic grating, Physica D 1
%   (1980), 267-290, especially Eqs. (21)-(22).
%   https://doi.org/10.1016/0167-2789(80)90026-3
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(s),s=1;end
validateattributes(s,{'numeric'}, ...
    {'scalar','integer','finite'});
if s<1,error('quantum_cat:Parameter','S must be a positive integer.');end
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('quantum_cat:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('quantum_cat',n);
    U.reference='Hannay-Berry, Physica D 1 (1980), 267-290, Eq. (21)';
    U.notes=['Quantization with periodic boundary conditions of ' ...
        'T=[2*S 1;4*S-1 2], S a positive integer. ' ...
        'Finite eigenangles depend arithmetically on N and S; ' ...
        'no generic eigenphase density is asserted.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('quantum_cat:Mode','Supported modes: matrix, operator, spectral.');
end
j=(0:n-1)';
qin=exp((2i*pi/n)*mod(s*mod(j.^2,n),n));
qout=exp((2i*pi/n)*mod(j.^2,n));
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_quantum_cat_apply(X,flag,n,qin,qout);
else
    F=exp((-2i*pi/n)*(j*j'))/sqrt(n);
    U=exp(1i*pi/4)*bsxfun(@times,qout, ...
        bsxfun(@times,F,qin.'));
end
end
