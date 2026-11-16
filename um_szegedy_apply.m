function Y=um_szegedy_apply(X,flag,P)
%UM_SZEGEDY_APPLY Matrix-free product of two Szegedy reflections.
%   State ordering is |i,j> at (i-1)*n+j; P is row-stochastic.
n=size(P,1);N=n*n;
if strcmpi(flag,'size'),Y=[N N];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_szegedy_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_szegedy_apply:Dimension','Input has wrong row count.');
end
a=sqrt(full(P));b=a.';
Y=zeros(size(X),'like',X);
for r=1:size(X,2)
    M=reshape(full(X(:,r)),n,n).'; % M(i,j) = X((i-1)*n+j).
    if adj
        M=reflect_b(M,b);
        M=reflect_a(M,a);
    else
        M=reflect_a(M,a);
        M=reflect_b(M,b);
    end
    Y(:,r)=reshape(M.',[],1);
end
end
function M=reflect_a(M,a)
v=sum(a.*M,2);
M=2*bsxfun(@times,a,v)-M;
end
function M=reflect_b(M,b)
v=sum(b.*M,1);
M=2*bsxfun(@times,b,v)-M;
end
