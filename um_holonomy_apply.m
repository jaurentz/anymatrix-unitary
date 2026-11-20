function Y=um_holonomy_apply(X,flag,P)
%UM_HOLONOMY_APPLY Sphere geodesic-loop rotations and adjoint.
n=size(P,1);m=size(P,2);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_holonomy_apply:Flag','Action flag must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_holonomy_apply:Dimension','X must have N rows.');
end
Y=X;
switch flag
    case 'notransp'
        for j=1:m
            Y=um_geodesic_apply(Y,'notransp',P(:,j),P(:,mod(j,m)+1),P(:,j)'*P(:,mod(j,m)+1));
        end
    case {'transp','ctransp'}
        for j=m:-1:1
            Y=um_geodesic_apply(Y,'transp',P(:,j),P(:,mod(j,m)+1),P(:,j)'*P(:,mod(j,m)+1));
        end
    otherwise,error('um_holonomy_apply:Flag','Unknown operator action.');
end
end
