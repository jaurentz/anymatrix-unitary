function Y=um_stiefel_apply(X,flag,frames,B)
%UM_STIEFEL_APPLY Efficient ambient Stiefel lift/adjoint action.
n=size(frames,1);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_stiefel_apply:Flag','Action flag must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n,n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_stiefel_apply:Dimension','X must have N rows.');
end
switch flag
    case 'notransp',C=B;
    case {'transp','ctransp'},C=B';
    otherwise,error('um_stiefel_apply:Flag','Unknown operator action.');
end
Y=X+frames*((C-eye(size(B,1)))*(frames'*X));
end
