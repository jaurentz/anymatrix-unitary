function Y=um_monomial_apply(X,flag,perm,weights)
%UM_MONOMIAL_APPLY Permutation and unimodular weights, including adjoint.
%   U(perm(j),j)=weights(j). No sparse or dense assembly is needed.
n=numel(perm);
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_monomial_apply:Dimension','Input must have n rows.');
end
switch lower(char(flag))
    case 'notransp'
        Y=X;
        Y(perm,:)=weights(:).*X;
    case {'transp','ctransp'}
        Y=conj(weights(:)).*X(perm,:);
    otherwise
        error('um_monomial_apply:Flag','Unknown operator action.');
end
end
