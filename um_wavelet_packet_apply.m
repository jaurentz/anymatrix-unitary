function Y=um_wavelet_packet_apply(X,flag,n,split,h,g)
%UM_WAVELET_PACKET_APPLY Orthogonal packet-tree analysis and adjoint.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_wavelet_packet_apply:Flag','Invalid operator flag.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_wavelet_packet_apply:Dimension','X must have N rows.');
end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_wavelet_packet_apply:Flag','Invalid operator action.');
end
if issparse(X),X=full(X);end
Y=X;nrhs=size(Y,2);
if strcmp(flag,'notransp')
    nodes=split(:)';
else
    nodes=fliplr(split(:)');
end
for node=nodes
    lev=floor(log2(node));
    m=n/2^lev;
    first=(node-2^lev)*m+1;
    rows=first+(0:m-1);
    half=m/2;
    if strcmp(flag,'notransp')
        A=zeros(half,nrhs,'like',Y);
        D=zeros(half,nrhs,'like',Y);
        for k=0:numel(h)-1
            idx=first+mod(2*(0:half-1)+k,m);
            A=A+h(k+1)*Y(idx,:);
            D=D+g(k+1)*Y(idx,:);
        end
        Y(rows,:)=[A;D];
    else
        A=Y(rows(1:half),:);
        D=Y(rows(half+1:end),:);
        restored=zeros(m,nrhs,'like',Y);
        for k=0:numel(h)-1
            local=mod(2*(0:half-1)+k,m)+1;
            restored(local,:)=restored(local,:)+h(k+1)*A+g(k+1)*D;
        end
        Y(rows,:)=restored;
    end
end
end
