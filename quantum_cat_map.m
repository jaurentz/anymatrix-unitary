function [U,properties]=quantum_cat_map(n,a,d,output)
%QUANTUM_CAT_MAP Finite metaplectic quantization of a hyperbolic torus map.
%   U=QUANTUM_CAT_MAP(N,A,D) quantizes the classical integral matrix
%       S = [A 1; A*D-1 D] in SL(2,Z),
%   with even positive A,D and A*D>=4 (hyperbolic, checkerboard parity).
%   Default A=D=2 gives [2 1;3 2]. For position labels 0<=j,k<N,
%      U(j+1,k+1) = exp(i*pi*(D*j^2-2*j*k+A*k^2)/N)/sqrt(N).
%   This has a specified global phase equal to 1; metaplectic lifts
%   are unique only up to a scalar. The phase periodicity is compatible
%   with integer indices for these even-diagonal classical maps, for
%   EVERY positive N. This is NOT a blanket formula for arbitrary
%   SL(2,Z) elements or the parity-incompatible [[2,1];[1,1]] map.
%   A diagonal chirp / DFT / diagonal chirp factorization makes U
%   exactly unitary. 'operator' uses O(N log N) FFT matvecs and adjoints.
%   No generic finite trace spectral density is implied by chaos.
%
%   J. H. Hannay and M. V. Berry, Physica D 1 (1980), 267-290.
%   https://doi.org/10.1016/0167-2789(80)90026-3
properties={'square','unitary','normal','complex','scalable', ...
            'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(a),a=2;end
if nargin<3||isempty(d),d=2;end
validateattributes(a,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if mod(a,2)~=0||mod(d,2)~=0||double(a)*double(d)<4
    error('quantum_cat_map:Classical', ...
          'A and D must be positive even integers with A*D>=4.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('quantum_cat_map:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('quantum_cat_map',n);
    U.reference='Hannay and Berry, Physica D 1 (1980), 267-290';
    U.notes=['An even-diagonal hyperbolic toral automorphism quantized ' ...
        'in a fixed projective phase convention. Its finite trace ' ...
        'eigenphase law is not presumed uniform.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('quantum_cat_map:Mode','Supported modes: matrix, operator, spectral.');
end
j=(0:n-1)';
left=exp(1i*pi*d*mod(j.^2,2*n)/n);
right=exp(1i*pi*a*mod(j.^2,2*n)/n);
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_cat_apply(X,flag,left,right);
    return
end
F=exp((-2i*pi/n)*(j*j'))/sqrt(n);
U=bsxfun(@times,left,bsxfun(@times,F,right.'));
end
