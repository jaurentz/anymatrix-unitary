function Y=um_wilson_apply(X,flag,links)
%UM_WILSON_APPLY Ordered product of unitary lattice-gauge link factors.
n=size(links,1);m=size(links,3);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_wilson_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_wilson_apply:Dimension','X must have N rows.');
end
Y=X;
switch flag
    case 'notransp'
        for j=m:-1:1,Y=links(:,:,j)*Y;end
    case {'transp','ctransp'}
        for j=1:m,Y=links(:,:,j)'*Y;end
    otherwise
        error('um_wilson_apply:Flag','Unknown action flag.');
end
end
