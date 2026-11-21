function Y=um_finite_koopman_apply(X,flag,perm,invperm)
%UM_FINITE_KOOPMAN_APPLY Exact permutation pullback / adjoint, O(N).
N=numel(perm);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_finite_koopman_apply:Flag','Invalid action flag.');
end
mode=lower(char(flag));
if strcmp(mode,'size'),Y=[N,N];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_finite_koopman_apply:Dimension','Input must have M^2 rows.');
end
switch mode
    case 'notransp',Y=X(perm,:);
    case {'transp','ctransp'},Y=X(invperm,:);
    otherwise
        error('um_finite_koopman_apply:Flag', ...
            'Action must be notransp, transp, ctransp, size.');
end
end
