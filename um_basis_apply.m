function Y=um_basis_apply(X,flag,V)
%UM_BASIS_APPLY Cached full spectral eigenbasis used by DPSS and graph FFT.
%   This is an operator-handle interface, NOT an asymptotically fast
%   action: V is a dense stored n-by-n orthonormal basis.
n=size(V,1);
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_basis_apply:Dimension','Input must have N rows.');
end
switch lower(char(flag))
    case 'notransp'
        Y=V'*X;
    case {'transp','ctransp'}
        Y=V*X;
    otherwise
        error('um_basis_apply:Flag','Unknown operator action.');
end
end
