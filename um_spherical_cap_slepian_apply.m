function Y=um_spherical_cap_slepian_apply(X,flag,C)
%UM_SPHERICAL_CAP_SLEPIAN_APPLY Real orthogonal block Slepian action.
%   C is computed by UM_SPHERICAL_CAP_COMMUTING, not by eig(K).
%   Each fixed-m Wigner-associated Legendre block diagonalizes the
%   exact cap concentration operator; only the SORTED permutations
%   couple the otherwise independent cos/sin block transforms.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_spherical_cap_slepian_apply:Flag','FLAG must be text.');
end
flag=lower(char(flag));
n=C.dimension;
if strcmp(flag,'size')
    Y=[n n];return;
end
if ~any(strcmp(flag,{'notransp','ctransp','transp'}))
    error('um_spherical_cap_slepian_apply:Flag', ...
        'Use notransp, ctransp, transp or size.');
end
if ~(isnumeric(X)&&ismatrix(X)&&size(X,1)==n)
    error('um_spherical_cap_slepian_apply:Dimension', ...
        'Input must have L^2 rows.');
end
Z=zeros(size(X),'like',X);
perm=C.sorted_permutation;
if strcmp(flag,'notransp')
    for j=1:numel(C.blocks)
        block=C.blocks{j};
        inds=block.ids;
        Z(inds,:)=block.V'*X(inds,:);
    end
    Y=Z(perm,:);
else
    Z(perm,:)=X;
    Y=zeros(size(X),'like',X);
    for j=1:numel(C.blocks)
        block=C.blocks{j};
        inds=block.ids;
        Y(inds,:)=block.V*Z(inds,:);
    end
end
end
