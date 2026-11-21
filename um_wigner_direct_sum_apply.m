function Y=um_wigner_direct_sum_apply(X,flag,blocks,N)
%UM_WIGNER_DIRECT_SUM_APPLY Degree-block SO(3) operator action.
%   BLOCKS{ell+1} is an explicit (2*ell+1)-square matrix OR a
%   handle F(X,flag) providing matrix-free actions. Numeric blocks
%   retain the original dense precomputation path. Handle blocks
%   allow sparse skew exponential-action with no O(L^4) setup.
%   For complex Wigner matrices, 'transp' means the NONCONJUGATE
%   transpose and 'ctransp' is the Hermitian adjoint; unlike the
%   old helper they are not conflated. Real blocks have both equal.
%   Numerical approximation/complexity follows each block's kernel.

if strcmpi(flag,'size'),Y=[N N];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_wigner_direct_sum_apply:Dimension', ...
        'Input must have N=L^2 rows.');
end
flag=lower(char(flag));
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_wigner_direct_sum_apply:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
Y=zeros(size(X),'like',X);
off=0;
for ell=0:numel(blocks)-1
    idx=off+(1:2*ell+1);
    block=blocks{ell+1};
    if isa(block,'function_handle')
        Y(idx,:)=block(X(idx,:),flag);
    elseif strcmp(flag,'notransp')
        Y(idx,:)=block*X(idx,:);
    elseif strcmp(flag,'transp')
        Y(idx,:)=block.'*X(idx,:);
    else
        Y(idx,:)=block'*X(idx,:);
    end
    off=off+2*ell+1;
end
end
