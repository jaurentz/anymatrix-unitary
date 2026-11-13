function Y=um_jaynes_apply(X,flag,a,b,d,edge1,edge2)
%UM_JAYNES_APPLY O(n) block action of the finite Jaynes-Cummings unitary.
%   a,b,d are n-vectors encoding the 2-by-2 blocks [a,b;b,d]
%   acting on rows (2,3),(4,5),..., (2*n,2*n+1).
%   The uncoupled rows 1 and 2*n+2 have phases edge1 and edge2.
n=numel(a);
N=2*(n+1);
if strcmpi(flag,'size'),Y=[N N];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_jaynes_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=N
    error('um_jaynes_apply:Dimension','Input must have 2*(n+1) rows.');
end
Y=complex(zeros(size(X)));
if adj
    Y(1,:)=conj(edge1)*X(1,:);
    Y(N,:)=conj(edge2)*X(N,:);
else
    Y(1,:)=edge1*X(1,:);
    Y(N,:)=edge2*X(N,:);
end
if n==0,return;end
i=(2:2:2*n)';
j=i+1;
if adj
    a=conj(a(:));b=conj(b(:));d=conj(d(:));
else
    a=a(:);b=b(:);d=d(:);
end
Y(i,:)=a.*X(i,:)+b.*X(j,:);
Y(j,:)=b.*X(i,:)+d.*X(j,:);
end
