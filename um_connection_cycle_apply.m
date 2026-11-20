function Y=um_connection_cycle_apply(X,flag,angles)
%UM_CONNECTION_CYCLE_APPLY O(N) rotation and permutation/adjoint.
N=numel(angles);n=2*N;
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_connection_cycle_apply:Flag','Action flag must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_connection_cycle_apply:Dimension','X must have 2N rows.');
end
c=cos(angles);s=sin(angles);
Y=zeros(size(X),'like',X);
j=(1:N)';d=mod(j,N)+1;
switch flag
    case 'notransp'
        Y(2*d-1,:)=c.*X(2*j-1,:)-s.*X(2*j,:);
        Y(2*d,:)=s.*X(2*j-1,:)+c.*X(2*j,:);
    case {'transp','ctransp'}
        Y(2*j-1,:)=c.*X(2*d-1,:)+s.*X(2*d,:);
        Y(2*j,:)= -s.*X(2*d-1,:)+c.*X(2*d,:);
    otherwise,error('um_connection_cycle_apply:Flag','Unknown operator action.');
end
end
