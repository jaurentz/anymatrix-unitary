function Y=um_butterfly_apply(X,flag,n,angles)
%UM_BUTTERFLY_APPLY O(N log N) binary-tree rotation butterfly action.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_butterfly_apply:Flag','Action flag must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_butterfly_apply:Dimension','Input must have N rows.');
end
switch flag
    case 'notransp'
        Y=apply_forward(X,angles,1);
    case {'transp','ctransp'}
        Y=apply_adjoint(X,angles,1);
    otherwise
        error('um_butterfly_apply:Flag','Unknown operator action.');
end
end

function Y=apply_forward(X,angles,k)
n=size(X,1);
if n==1,Y=X;return;end
m=n/2;c=cos(angles(k));s=sin(angles(k));
L=apply_forward(X(1:m,:),angles,k+1);
R=apply_forward(X(m+1:end,:),angles,k+m);
Y=[c*L+s*R;-s*L+c*R];
end

function Y=apply_adjoint(X,angles,k)
n=size(X,1);
if n==1,Y=X;return;end
m=n/2;c=cos(angles(k));s=sin(angles(k));
L=apply_adjoint(c*X(1:m,:)-s*X(m+1:end,:),angles,k+1);
R=apply_adjoint(s*X(1:m,:)+c*X(m+1:end,:),angles,k+m);
Y=[L;R];
end
