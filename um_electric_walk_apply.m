function Y=um_electric_walk_apply(X,flag,d,L,M)
%UM_ELECTRIC_WALK_APPLY Sparse O(N) phase/CMV-factor quantum walk.
%   W=D*L*M. Supports distinct ordinary transpose and adjoint
%   for the complex electric field and twisted boundaries.
N=numel(d);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_electric_walk_apply:Flag','Action must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[N,N];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_electric_walk_apply:Dimension', ...
       'X must have 2N rows.');
end
switch flag
    case 'notransp'
        Y=bsxfun(@times,d,L*(M*X));
    case 'ctransp'
        Y=M'*(L'*bsxfun(@times,conj(d),X));
    case 'transp'
        Y=M.'*(L.'*bsxfun(@times,d,X));
    otherwise
        error('um_electric_walk_apply:Flag', ...
            'Action must be notransp, transp, ctransp or size.');
end
end
