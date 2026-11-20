function Y=um_hyperbolic_holonomy_apply(X,flag,n,theta)
%UM_HYPERBOLIC_HOLONOMY_APPLY Two-dimensional holonomy plus spectators.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_hyperbolic_holonomy_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_hyperbolic_holonomy_apply:Dimension','X must have N rows.');
end
switch flag
    case 'notransp'
        s=1;
    case {'transp','ctransp'}
        s=-1;
    otherwise
        error('um_hyperbolic_holonomy_apply:Flag','Unknown action flag.');
end
c=cos(theta);z=sin(theta);
Y=X;
Y(1,:)=c*X(1,:)-s*z*X(2,:);
Y(2,:)=s*z*X(1,:)+c*X(2,:);
end
