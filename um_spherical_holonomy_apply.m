function Y=um_spherical_holonomy_apply(X,flag,vertices,c)
%UM_SPHERICAL_HOLONOMY_APPLY Factorized closed geodesic polygon transport.
n=size(vertices,1);m=size(vertices,2);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_spherical_holonomy_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_spherical_holonomy_apply:Dimension','X must have N rows.');
end
Y=X;
switch flag
    case 'notransp'
        for j=1:m
            y=vertices(:,mod(j,m)+1);
            Y=um_geodesic_apply(Y,'notransp',vertices(:,j),y,c(j));
        end
    case {'transp','ctransp'}
        for j=m:-1:1
            y=vertices(:,mod(j,m)+1);
            Y=um_geodesic_apply(Y,'transp',vertices(:,j),y,c(j));
        end
    otherwise
        error('um_spherical_holonomy_apply:Flag','Unknown action flag.');
end
end
