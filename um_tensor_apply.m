function Y=um_tensor_apply(X,flag,factors)
%UM_TENSOR_APPLY Matrix-free action of ordered Kronecker factors.
%   FACTORS={A1,...,At} represents kron(...kron(A1,A2),...,At).
%   The modewise algorithm uses no full Kronecker matrix and accepts
%   multiple right-hand sides. 'transp' denotes conjugate transpose.
t=numel(factors);
dims=zeros(1,t);
for j=1:t
    A=factors{j};
    if ~isnumeric(A)||ndims(A)>2||size(A,1)~=size(A,2)
        error('um_tensor_apply:Factor','All factors must be square numeric matrices.');
    end
    dims(j)=size(A,1);
end
n=prod(dims);
if strcmpi(flag,'size'),Y=[n n];return;end
if any(strcmpi(flag,{'transp','ctransp'}))
    adj=true;
elseif strcmpi(flag,'notransp')
    adj=false;
else
    error('um_tensor_apply:Flag','Invalid operator action.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_tensor_apply:Dimension','Input must have prod(dims) rows.');
end
Y=full(X);
nrhs=size(Y,2);
for j=1:t
    right=prod(dims(j+1:end));
    left=prod(dims(1:j-1));
    B=reshape(Y,[right,dims(j),left*nrhs]);
    B=permute(B,[2 1 3]);
    B=reshape(B,dims(j),[]);
    if adj
        B=factors{j}'*B;
    else
        B=factors{j}*B;
    end
    B=reshape(B,[dims(j),right,left*nrhs]);
    B=ipermute(B,[2 1 3]);
    Y=reshape(B,n,nrhs);
end
end
