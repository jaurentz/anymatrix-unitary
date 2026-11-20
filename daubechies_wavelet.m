function [U,properties]=daubechies_wavelet(n,levels,output,order)
%DAUBECHIES_WAVELET Periodic finite orthonormal Daubechies analysis transform.
%   U=DAUBECHIES_WAVELET(N,LEVELS,OUTPUT,ORDER) builds an N-by-N real
%   orthogonal multilevel wavelet analysis matrix with periodic boundary
%   handling. N must be a power of two (including N=1).
%   ORDER=2,3,4 selects db2, db3, db4, respectively (filter lengths
%   4,6,8); default ORDER=2. ORDER is the optional FOURTH argument to
%   preserve the legacy DAUBECHIES_WAVELET(N,LEVELS,OUTPUT) interface.
%   LEVELS defaults to log2(N); LEVELS=0
%   returns the identity.
%
%   At a level of current size m, the row-index j=0:m/2-1 outputs are
%      A(j)=sum_k h(k)*x(mod(2*j+k,m)),
%      D(j)=sum_k g(k)*x(mod(2*j+k,m)),
%   with k zero-based, h the standard extremal-phase reconstruction
%   low-pass filter, and g(k)=(-1)^k*h(L-1-k). The signs and periodization
%   therefore differ by shifts/flips from some DWT package conventions.
%   The stage is orthogonal; successive stages act on approximation
%   coefficients only, so U is exactly orthogonal in exact arithmetic.
%   Explicit mode materializes U; 'operator' performs O(N*ORDER*LEVELS)
%   forward/adjoint actions without building a matrix.
%
%   [15] I. Daubechies, Ten Lectures on Wavelets, SIAM (1992).
%   https://doi.org/10.1137/1.9781611970104
%   Tabulated coefficients independently checked against
%   https://github.com/PyWavelets/pywt/blob/main/pywt/_extensions/c/wavelets_coeffs.template.h
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
maxlevels=round(log2(n));
if 2^maxlevels~=n
    error('daubechies_wavelet:Order','N must be a power of two.');
end
if nargin<2||isempty(levels),levels=maxlevels;end
validateattributes(levels,{'numeric'}, ...
    {'scalar','integer','finite'});
if levels<0 || levels>maxlevels
    error('daubechies_wavelet:Levels','LEVELS must be between 0 and log2(N).');
end
if nargin<3||isempty(output),output='matrix';end
if nargin<4||isempty(order),order=2;end
validateattributes(order,{'numeric'}, ...
    {'scalar','integer','finite'});
if order<2 || order>4
    error('daubechies_wavelet:Filter','ORDER must be 2, 3, or 4.');
end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('daubechies_wavelet:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('daubechies_wavelet',n);
    U.reference='Daubechies, Ten Lectures on Wavelets (1992)';
    if levels==0
        U.finite.kind='exact';
        U.finite.eigenvalues=@()ones(n,1);
        U.finite.moments=@(k)ones(size(k));
        U.notes='Zero decomposition levels: exactly the identity matrix.';
    else
        U.notes=['Orthogonal finite periodic wavelet analysis, order db2, db3, ' ...
            'or db4. The finite trace eigenphase law is not the scaling ' ...
            'filter frequency response and has no generic closed form.'];
    end
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('daubechies_wavelet:Mode', ...
        'Supported modes: matrix, operator, spectral.');
end
h=db_lowpass(order);
g=(-1).^(0:numel(h)-1)'.*flipud(h);
F=@(X,flag)um_daubechies_apply(X,flag,n,levels,h,g);
if any(strcmp(output,{'operator','implicit'}))
    U=F;
else
    U=F(eye(n),'notransp');
end
end

function h=db_lowpass(order)
% Extremal-phase reconstruction low-pass coefficients from PyWavelets.
switch order
    case 2
        h=[ ...
          .4829629131445341433748715998644486838; ...
          .8365163037378079055752937809168732035; ...
          .2241438680420133810259727622404003555; ...
         -.1294095225512603811744494188120241642];
    case 3
        h=[ ...
          .3326705529500826159985115891390056300; ...
          .8068915093110925764944936040887134905; ...
          .4598775021184915700951519421476167208; ...
         -.1350110200102545886963899066993744806; ...
         -.0854412738820266616928191691817733115; ...
          .0352262918857095366027406647155100293];
    case 4
        h=[ ...
          .2303778133088965008632911830440708500; ...
          .7148465705529156470899219552739926037; ...
          .6308807679298589078817163383006152202; ...
         -.0279837694168598542114137471800753854; ...
         -.1870348117190930840795706727890814196; ...
          .0308413818355607636272193625349590502; ...
          .0328830116668851997354075135492443887; ...
         -.0105974017850690321048832085240272207];
end
end
