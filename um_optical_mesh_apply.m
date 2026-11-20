function Y=um_optical_mesh_apply(X,flag,n,pairs,c,s,phase,ext)
%UM_OPTICAL_MESH_APPLY Direct core-mesh forward and Hermitian adjoint.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_optical_mesh_apply:Flag','Invalid flag.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_optical_mesh_apply:Dimension','Expected N input rows.');
end
if issparse(X),X=full(X);end
Y=complex(X);
switch flag
    case 'notransp'
        for t=1:numel(pairs)
            j=pairs(t);a=Y(j,:);b=Y(j+1,:);
            Y(j,:)=phase(t)*c(t)*a-s(t)*b;
            Y(j+1,:)=phase(t)*s(t)*a+c(t)*b;
        end
        Y=bsxfun(@times,ext,Y);
    case {'transp','ctransp'}
        Y=bsxfun(@times,conj(ext),Y);
        for t=numel(pairs):-1:1
            j=pairs(t);a=Y(j,:);b=Y(j+1,:);
            Y(j,:)=conj(phase(t))*(c(t)*a+s(t)*b);
            Y(j+1,:)=-s(t)*a+c(t)*b;
        end
    otherwise
        error('um_optical_mesh_apply:Flag','Unknown operator flag.');
end
end
