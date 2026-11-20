function Y=um_quantum_cat_apply(X,flag,n,qin,qout)
%UM_QUANTUM_CAT_APPLY FFT of chirp factors for Hannay-Berry propagator.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_quantum_cat_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_quantum_cat_apply:Dimension','Input must have N rows.');
end
switch flag
    case 'notransp'
        Y=exp(1i*pi/4)*bsxfun(@times,qout, ...
            fft(bsxfun(@times,qin,X),[],1)/sqrt(n));
    case {'transp','ctransp'}
        Y=exp(-1i*pi/4)*bsxfun(@times,conj(qin), ...
            ifft(bsxfun(@times,conj(qout),X),[],1)*sqrt(n));
    otherwise
        error('um_quantum_cat_apply:Flag','Unknown operator action.');
end
end
