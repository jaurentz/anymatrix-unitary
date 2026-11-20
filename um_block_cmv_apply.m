function Y=um_block_cmv_apply(X,flag,n,q,cores,beta)
%UM_BLOCK_CMV_APPLY Local-factor block-CMV forward/adjoint.
N=n*q;
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_block_cmv_apply:Flag','FLAG must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[N N];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_block_cmv_apply:Dimension','Expected N*Q input rows.');
end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_block_cmv_apply:Flag','Unknown action.');
end
if issparse(X),X=full(X);end
Y=X;
last=(n-1)*q+(1:q);
if strcmp(flag,'notransp')
    % U=L*M: apply even M cores, then odd L cores.
    for k=2:2:n-1
        rows=(k-1)*q+(1:2*q);
        Y(rows,:)=cores{k}*Y(rows,:);
    end
    if mod(n,2)==0,Y(last,:)=beta'*Y(last,:);end
    for k=1:2:n-1
        rows=(k-1)*q+(1:2*q);
        Y(rows,:)=cores{k}*Y(rows,:);
    end
    if mod(n,2)==1,Y(last,:)=beta'*Y(last,:);end
else
    % U'=M'*L', with reversal of factor order and core adjoints.
    if mod(n,2)==1,Y(last,:)=beta*Y(last,:);end
    for k=1:2:n-1
        rows=(k-1)*q+(1:2*q);
        Y(rows,:)=cores{k}'*Y(rows,:);
    end
    if mod(n,2)==0,Y(last,:)=beta*Y(last,:);end
    for k=2:2:n-1
        rows=(k-1)*q+(1:2*q);
        Y(rows,:)=cores{k}'*Y(rows,:);
    end
end
end
