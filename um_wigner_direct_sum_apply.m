function Y=um_wigner_direct_sum_apply(X,flag,blocks,N)
%UM_WIGNER_DIRECT_SUM_APPLY Precomputed irreducible SO(3) block action.
%   BLOCKS{ell+1} is (2*ell+1)-square unitary/orthogonal.
%   For L blocks, N=L^2, storage and multiply O(sum (2ell+1)^2)
%   = O(L^3) = O(N^(3/2)). Precomputing blocks uses dense EXPM,
%   O(sum ell^3)=O(L^4), so no O(L^2 log L) "fast Wigner"
%   or Gimbutas--Greengard pseudospectral algorithm is claimed.
if strcmpi(flag,'size'),Y=[N N];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_wigner_direct_sum_apply:Dimension', ...
        'Input must have N=L^2 rows.');
end
switch lower(char(flag))
    case 'notransp'
        adj=false;
    case {'transp','ctransp'}
        adj=true;
    otherwise
        error('um_wigner_direct_sum_apply:Flag', ...
            'Use notransp, transp, ctransp or size.');
end
Y=zeros(size(X),'like',X);
off=0;
for ell=0:numel(blocks)-1
    idx=off+(1:2*ell+1);
    if adj,Y(idx,:)=blocks{ell+1}'*X(idx,:);
    else,Y(idx,:)=blocks{ell+1}*X(idx,:);end
    off=off+2*ell+1;
end
end
