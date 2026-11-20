function Y=um_cat_apply(X,flag,left,right)
%UM_CAT_APPLY Chirp-DFT-chirp action for parity-admissible cat maps.
n=numel(left);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_cat_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_cat_apply:Dimension','X must have N rows.');
end
switch flag
    case 'notransp'
        Y=bsxfun(@times,right,X);
        Y=fft(Y,[],1)/sqrt(n);
        Y=bsxfun(@times,left,Y);
    case {'transp','ctransp'}
        Y=bsxfun(@times,conj(left),X);
        Y=ifft(Y,[],1)*sqrt(n);
        Y=bsxfun(@times,conj(right),Y);
    otherwise
        error('um_cat_apply:Flag','Unknown action.');
end
end
