function Y=um_dirac_apply(X,flag,n,nu,m)
%UM_DIRAC_APPLY Matrix-free forward and conjugate-transpose Dirac actions.
%   Inputs are ordered [up(1:n,:); down(1:n,:)].
if strcmpi(flag,'size'),Y=[2*n 2*n];return;end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=2*n
    error('um_dirac_apply:Dimension','Input must have 2*n rows.');
end
up=X(1:n,:);down=X(n+1:end,:);
right=[n 1:n-1];left=[2:n 1];
switch lower(char(flag))
    case 'notransp'
        Y=[nu*up(right,:)-1i*m*down; ...
           -1i*m*up+nu*down(left,:)];
    case {'transp','ctransp'}
        Y=[nu*up(left,:)+1i*m*down; ...
           1i*m*up+nu*down(right,:)];
    otherwise
        error('um_dirac_apply:Flag', ...
            'Use notransp, transp, ctransp or size.');
end
end
