function Y=um_affine_apply(X,flag,d1,d2)
%UM_AFFINE_APPLY Chirp--DFT--chirp using the shared Fourier kernel.
%   A=diag(d2)*F*diag(d1), with unit-modulus chirp factors.
n=numel(d1);
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_affine_apply:Dimension','Input must have n rows.');
end
switch lower(char(flag))
    case 'notransp'
        Y=d2(:).*um_fourier_apply(d1(:).*X,'notransp',0);
    case {'transp','ctransp'}
        Y=conj(d1(:)).*um_fourier_apply(conj(d2(:)).*X,'transp',0);
    otherwise
        error('um_affine_apply:Flag','Unknown transform action.');
end
end
