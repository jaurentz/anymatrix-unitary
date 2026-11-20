function Y=um_periodic_sinc_apply(X,flag,phase,n)
%UM_PERIODIC_SINC_APPLY FFT-based periodic fractional shift and adjoint.
if isstring(flag),flag=char(flag);end
if ~ischar(flag)
    error('um_periodic_sinc_apply:Flag','Invalid operator flag.');
end
if strcmpi(flag,'size')
    Y=[n n];return
end
if size(X,1)~=n
    error('um_periodic_sinc_apply:Dimension','Expected N input rows.');
end
switch lower(flag)
    case 'notransp'
        d=phase;
    case {'transp','ctransp'}
        d=conj(phase);
    otherwise
        error('um_periodic_sinc_apply:Flag','Invalid operator flag.');
end
Y=ifft(bsxfun(@times,d,fft(X,[],1)),[],1);
if isreal(X),Y=real(Y);end
end
