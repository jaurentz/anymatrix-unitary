function Y=um_gate_apply(X,gate,sites,n)
%UM_GATE_APPLY Apply a 2-by-2 or 4-by-4 gate on an n-qubit state array.
%   Site 1 is the MOST-significant bit, as in QM_TWO_SITE_GATE.
%   SITES order is significant for a two-site gate, including wrap bonds.
k=numel(sites);
if k<1 || k>2 || numel(unique(sites))~=k || any(sites<1) || any(sites>n)
    error('um_gate_apply:Sites','Sites must be one or two distinct qubits.');
end
if ~isequal(size(gate),[2^k 2^k])
    error('um_gate_apply:Gate','Gate dimension does not match sites.');
end
if size(X,1)~=2^n,error('um_gate_apply:Dimension','Invalid state dimension.');end
nc=size(X,2);
axes=n-fliplr(sites)+1;  % MATLAB first tensor axis is least significant.
p=[axes,setdiff(1:n,axes,'stable'),n+1];
T=permute(reshape(X,[2*ones(1,n),nc]),p);
T=gate*reshape(T,2^k,[]);
T=reshape(T,[2*ones(1,n),nc]);
Y=reshape(ipermute(T,p),2^n,nc);
end
