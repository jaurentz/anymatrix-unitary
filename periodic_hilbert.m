function [U,properties] = periodic_hilbert(n,dc,nyquist)
%PERIODIC_HILBERT Real orthogonal completion of the periodic Hilbert filter.
%   The Fourier multiplier is -1i at positive frequencies, +1i at
%   negative frequencies, with independently assigned DC and (for even N)
%   Nyquist signs +1/-1. Defaults DC=NYQUIST=+1.
%   Without these unimodular endpoint signs the conventional discrete
%   Hilbert transform is not unitary on the entire finite signal space.
%   Eigenvalues are exactly the specified Fourier multipliers.
properties={'square','unitary','normal','orthogonal','real', ...
    'circulant','scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(dc),dc=1;end
if nargin<3 || isempty(nyquist),nyquist=1;end
if ~(isnumeric(dc)&&isscalar(dc)&&isreal(dc)&&any(dc==[-1,1])) || ...
   ~(isnumeric(nyquist)&&isscalar(nyquist)&&isreal(nyquist)&&any(nyquist==[-1,1]))
    error('periodic_hilbert:Signs','DC and NYQUIST must be +1 or -1.');
end
d=zeros(n,1);d(1)=dc;
d(2:ceil(n/2))=-1i;
d(floor(n/2)+2:end)=1i;
if mod(n,2)==0,d(n/2+1)=nyquist;end
column=ifft(d);
idx=mod((0:n-1)'-(0:n-1),n)+1;
U=real(column(idx));
end
