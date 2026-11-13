function Y=um_zak_apply(X,flag,p,q)
%UM_ZAK_APPLY Finite Zak FFTs on length-q polyphase columns.
%   Z=kron(F_q,I_p); the index j runs fastest in j+p*k.
n=p*q;
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_zak_apply:Dimension','Input must have n rows.');
end
m=size(X,2);
V=reshape(X,p,q,m);
switch lower(char(flag))
    case 'notransp'
        V=fft(V,[],2)/sqrt(q);
    case {'transp','ctransp'}
        V=ifft(V,[],2)*sqrt(q);
    otherwise
        error('um_zak_apply:Flag','Unknown transform action.');
end
Y=reshape(V,n,m);
end
