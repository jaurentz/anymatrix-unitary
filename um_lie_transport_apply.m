function Y=um_lie_transport_apply(X,flag,n,alpha)
%UM_LIE_TRANSPORT_APPLY Orthogonal pair rotations in Lie-algebra coordinates.
d=n*n;
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_lie_transport_apply:Flag','Action flag must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[d d];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=d
    error('um_lie_transport_apply:Dimension','X must have N^2 rows.');
end
switch flag
    case 'notransp',sgn=1;
    case {'transp','ctransp'},sgn=-1;
    otherwise,error('um_lie_transport_apply:Flag','Unknown operator action.');
end
a=n+2*(0:numel(alpha)-1)'+1;b=a+1;
c=cos(alpha);s=sin(alpha)*sgn;
Y=X;
Y(a,:)=c.*X(a,:)+s.*X(b,:);
Y(b,:)=-s.*X(a,:)+c.*X(b,:);
end
