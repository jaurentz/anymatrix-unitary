function Y=um_daubechies_apply(X,flag,n,levels,h,g)
%UM_DAUBECHIES_APPLY Forward and adjoint periodic wavelet filter bank.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_daubechies_apply:Flag','FLAG must be a scalar string.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_daubechies_apply:Dimension','Input must have N rows.');
end
Y=X;nrhs=size(X,2);
switch flag
    case 'notransp'
        for lev=1:levels
            m=n/2^(lev-1);
            half=m/2;
            A=zeros(half,nrhs,'like',Y);
            D=zeros(half,nrhs,'like',Y);
            for k=0:numel(h)-1
                idx=mod(2*(0:half-1)+k,m)+1;
                A=A+h(k+1)*Y(idx,:);
                D=D+g(k+1)*Y(idx,:);
            end
            Y(1:m,:)=[A;D];
        end
    case {'transp','ctransp'}
        for lev=levels:-1:1
            m=n/2^(lev-1);
            half=m/2;
            A=Y(1:half,:);
            D=Y(half+1:m,:);
            restored=zeros(m,nrhs,'like',Y);
            for k=0:numel(h)-1
                idx=mod(2*(0:half-1)+k,m)+1;
                restored(idx,:)=restored(idx,:)+h(k+1)*A+g(k+1)*D;
            end
            Y(1:m,:)=restored;
        end
    otherwise
        error('um_daubechies_apply:Flag','Unknown operator action.');
end
end
