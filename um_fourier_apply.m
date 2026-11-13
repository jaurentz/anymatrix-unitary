function Y=um_fourier_apply(X,flag,shift)
%UM_FOURIER_APPLY Normalized shifted Fourier transform and its adjoint.
%   F(j,k)=exp(-2*pi*1i*(j+shift)*(k+shift)/n)/sqrt(n).
if nargin<3,shift=0;end
n=size(X,1);q=(0:n-1)';
d=exp(-2*pi*1i*shift*q/n);
c=exp(-2*pi*1i*shift^2/n);
switch lower(char(flag))
    case 'notransp'
        Y=(c/sqrt(n))*(d.*fft(d.*X,[],1));
    case {'transp','ctransp'}
        Y=(conj(c)*sqrt(n))*(conj(d).*ifft(conj(d).*X,[],1));
    otherwise
        error('um_fourier_apply:Flag','Unknown Fourier action.');
end
end
