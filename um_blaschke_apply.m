function Y=um_blaschke_apply(X,flag,n,b,V,U0)
%UM_BLASCHKE_APPLY Matrix-free application of rank-one Potapov factors.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_blaschke_apply:Flag','Action flag must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_blaschke_apply:Dimension','Input must have N rows.');
end
switch flag
    case 'notransp'
        Y=X;
        for j=numel(b):-1:1
            v=V(:,j);
            Y=Y+(b(j)-1)*v*(v'*Y);
        end
        Y=U0*Y;
    case {'transp','ctransp'}
        Y=U0'*X;
        for j=1:numel(b)
            v=V(:,j);
            Y=Y+(conj(b(j))-1)*v*(v'*Y);
        end
    otherwise
        error('um_blaschke_apply:Flag','Unknown operator action.');
end
end
