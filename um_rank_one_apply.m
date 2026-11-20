function Y=um_rank_one_apply(X,flag,H,n,v,tau)
%UM_RANK_ONE_APPLY Core-product action composed with a unitary reflector.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_rank_one_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_rank_one_apply:Dimension','X must have N rows.');
end
switch flag
    case 'notransp'
        Z=X+tau*v*(v'*X);
        Y=H(Z,'notransp');
    case {'transp','ctransp'}
        Z=H(X,'transp');
        Y=Z+conj(tau)*v*(v'*Z);
    otherwise
        error('um_rank_one_apply:Flag','Unknown action flag.');
end
end
