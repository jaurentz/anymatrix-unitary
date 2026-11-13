function Y=um_dft_apply(X,flag,n)
%UM_DFT_APPLY Expose the shared unitary Fourier action with size metadata.
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_dft_apply:Dimension','Input must have N rows.');
end
switch lower(char(flag))
    case 'notransp'
        Y=um_fourier_apply(X,'notransp',0);
    case {'transp','ctransp'}
        Y=um_fourier_apply(X,'transp',0);
    otherwise
        error('um_dft_apply:Flag','Unknown operator action.');
end
end
