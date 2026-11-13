function Y=um_fft_apply(X,flag,kind,a,b)
%UM_FFT_APPLY FFT-based circulant, kicked and baker unitary actions.
%   circulant: U=F'*diag(a)*F
%   kick: U=diag(b)*F'*diag(a)*F
%   baker: a=n, b=shift, U=F_n'*blkdiag(F_(n/2),F_(n/2))
switch lower(char(kind))
    case {'circulant','kick'},n=numel(a);
    case 'baker',n=a;
    otherwise,error('um_fft_apply:Kind','Unknown FFT operator kind.');
end
if strcmpi(flag,'size'),Y=[n n];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_fft_apply:Flag','Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_fft_apply:Dimension','Input must have n rows.');
end
% MATLAB FFT routines do not accept sparse arrays. The FFT result is
% generally dense, so expand sparse input only after checking dimensions.
% This also permits canonical-basis probes used for exact trace estimates.
if issparse(X),X=full(X);end
switch lower(char(kind))
    case 'circulant'
        if adj,phases=conj(a(:));else,phases=a(:);end
        Y=ifft(phases.*fft(X,[],1),[],1);
    case 'kick'
        if adj
            Y=ifft(conj(a(:)).*fft(conj(b(:)).*X,[],1),[],1);
        else
            Y=b(:).*ifft(a(:).*fft(X,[],1),[],1);
        end
    case 'baker'
        h=n/2;shift=b;
        if adj
            T=um_fourier_apply(X,'notransp',shift);
            Y=[um_fourier_apply(T(1:h,:),'transp',shift); ...
               um_fourier_apply(T(h+1:end,:),'transp',shift)];
        else
            T=[um_fourier_apply(X(1:h,:),'notransp',shift); ...
               um_fourier_apply(X(h+1:end,:),'notransp',shift)];
            Y=um_fourier_apply(T,'transp',shift);
        end
end
end
