function Y=um_walsh_apply(X,flag,perm)
%UM_WALSH_APPLY Fast normalized Sylvester butterfly with row ordering.
%   U=P*H, where rows perm = sequency or natural ordering.
n=numel(perm);
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_walsh_apply:Dimension','Input must have N rows.');
end
switch lower(char(flag))
    case 'notransp'
        Y=X;
        adj=false;
    case {'transp','ctransp'}
        Y=zeros(size(X),'like',X);
        Y(perm,:)=X;
        adj=true;
    otherwise
        error('um_walsh_apply:Flag','Unknown operator action.');
end
m=size(Y,2);
h=1;
while h<n
    B=reshape(Y,2*h,n/(2*h),m);
    a=B(1:h,:,:);b=B(h+1:2*h,:,:);
    Y=reshape(cat(1,a+b,a-b),n,m);
    h=2*h;
end
Y=Y/sqrt(n);
if ~adj,Y=Y(perm,:);end
end
