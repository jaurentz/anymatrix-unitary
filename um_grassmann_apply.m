function Y=um_grassmann_apply(X,flag,B,angles)
%UM_GRASSMANN_APPLY Compact O(Nr) unitary geodesic lift.
n=size(B,1);r=numel(angles);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_grassmann_apply:Flag','Action flag must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_grassmann_apply:Dimension','X must have N rows.');
end
switch flag
    case 'notransp',sgn=1;
    case {'transp','ctransp'},sgn=-1;
    otherwise,error('um_grassmann_apply:Flag','Unknown operator action.');
end
T=B'*X;C=cos(angles)-1;S=sin(angles)*sgn;
Y=X+B*[C.*T(1:r,:)-S.*T(r+1:end,:); ...
        S.*T(1:r,:)+C.*T(r+1:end,:)];
end
