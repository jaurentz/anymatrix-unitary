function Y=um_bundle_transport_apply(X,flag,permutation,links)
%UM_BUNDLE_TRANSPORT_APPLY Fast block permutation/adjoint.
n=numel(permutation);d=size(links,1);N=n*d;
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_bundle_transport_apply:Flag','Invalid flag.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[N N];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_bundle_transport_apply:Dimension','Expected N*d input rows.');
end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_bundle_transport_apply:Flag','Invalid flag.');
end
if issparse(X),X=full(X);end
Y=zeros(N,size(X,2));
if ~isreal(X)||~isreal(links),Y=complex(Y);end
for j=1:n
    from=(j-1)*d+(1:d);
    to=(permutation(j)-1)*d+(1:d);
    if strcmp(flag,'notransp')
        Y(to,:)=links(:,:,j)*X(from,:);
    else
        Y(from,:)=links(:,:,j)'*X(to,:);
    end
end
end
