function Y=um_sphere_apply(X,flag,p,q)
%UM_SPHERE_APPLY O(N) exact-rank-two sphere-rotation/adjoint action.
n=numel(p);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_sphere_apply:Flag','Action flag must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_sphere_apply:Dimension','X must have N rows.');
end
switch flag
    case 'notransp',sgn=1;
    case {'transp','ctransp'},sgn=-1;
    otherwise,error('um_sphere_apply:Flag','Unknown operator action.');
end
A=q*(p'*X)-p*(q'*X);
A2=q*(p'*A)-p*(q'*A);
Y=X+sgn*A+A2/(1+p'*q);
end
