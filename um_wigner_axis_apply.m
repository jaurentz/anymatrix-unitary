function Y=um_wigner_axis_apply(X,flag,phase,perm)
%UM_WIGNER_AXIS_APPLY Exact monomial Wigner z/y-pi Euler rotation.
%   U(row,perm(row))=phase(row).  Row-wise phases are unimodular.
%   Unlike a generic dense Wigner d, both U and U' are O(n) actions.
n=numel(phase);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_wigner_axis_apply:Flag','FLAG must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_wigner_axis_apply:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_wigner_axis_apply:Dimension', ...
        'Input must have 2*ell+1 rows.');
end
if strcmp(flag,'notransp')
    Y=phase(:).*X(perm,:);
else
    Y=zeros(size(X),'like',X);
    if strcmp(flag,'ctransp'),p=conj(phase(:));
    else,p=phase(:);end
    Y(perm,:)=p.*X;
end
end
