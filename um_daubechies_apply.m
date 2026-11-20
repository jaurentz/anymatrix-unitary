function Y=um_daubechies_apply(X,flag,n,levels)
%UM_DAUBECHIES_APPLY Periodized D4 2-channel filter, forward and adjoint.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_daubechies_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_daubechies_apply:Dimension','Input must have N rows.');
end
h=[1+sqrt(3),3+sqrt(3),3-sqrt(3),1-sqrt(3)]/(4*sqrt(2));
g=[h(4),-h(3),h(2),-h(1)];
Y=X;
switch flag
    case 'notransp'
        for ell=1:levels
            m=n/2^(ell-1);
            Y(1:m,:)=stage(Y(1:m,:),m,h,g,false);
        end
    case {'transp','ctransp'}
        for ell=levels:-1:1
            m=n/2^(ell-1);
            Y(1:m,:)=stage(Y(1:m,:),m,h,g,true);
        end
    otherwise
        error('um_daubechies_apply:Flag','Unknown action.');
end
end

function Z=stage(X,m,h,g,adjoint)
nrhs=size(X,2); half=m/2;
if ~adjoint
    Z=zeros(m,nrhs,'like',X);
    for j=0:half-1
        idx=mod(2*j+(0:3),m)+1;
        Z(j+1,:)=h*X(idx,:);
        Z(half+j+1,:)=g*X(idx,:);
    end
else
    Z=zeros(m,nrhs,'like',X);
    for j=0:half-1
        idx=mod(2*j+(0:3),m)+1;
        for k=1:4
            Z(idx(k),:)=Z(idx(k),:)+h(k)*X(j+1,:) ...
                         +g(k)*X(half+j+1,:);
        end
    end
end
end
