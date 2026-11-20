function Y=um_hadamard_apply(X,flag,n,m,twists)
%UM_HADAMARD_APPLY Fast two-axis Dita-Fourier Hadamard matrix actions.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_hadamard_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_hadamard_apply:Dimension','X must have N rows.');
end
r=n/m;nrhs=size(X,2);
Z=reshape(X,r,m,nrhs);
P=exp(1i*twists);
switch flag
    case 'notransp'
        Z=fft(Z,[],1)/sqrt(r);
        Z=bsxfun(@times,Z,P);
        Z=fft(Z,[],2)/sqrt(m);
    case {'transp','ctransp'}
        Z=ifft(Z,[],2)*sqrt(m);
        Z=bsxfun(@times,Z,conj(P));
        Z=ifft(Z,[],1)*sqrt(r);
    otherwise
        error('um_hadamard_apply:Flag','Unknown action flag.');
end
Y=reshape(Z,n,nrhs);
end
