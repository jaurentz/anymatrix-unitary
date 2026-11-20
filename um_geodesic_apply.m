function Y=um_geodesic_apply(X,flag,x,y,c)
%UM_GEODESIC_APPLY Matrix-free forward/adjoint minimal SO(n) rotation.
n=numel(x);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_geodesic_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_geodesic_apply:Dimension','Input must have N rows.');
end
switch flag
    case 'notransp'
        s=1;
    case {'transp','ctransp'}
        s=-1;
    otherwise
        error('um_geodesic_apply:Flag','Unknown action flag.');
end
ax=x'*X; ay=y'*X;
Y=X+s*(y*ax-x*ay) ...
    +(c*(y*ax+x*ay)-y*ay-x*ax)/(1+c);
end
