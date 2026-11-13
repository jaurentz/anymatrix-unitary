function Y=um_circuit_apply(X,flag,n,gates,depth)
%UM_CIRCUIT_APPLY Chronological local-gate circuit and conjugate adjoint.
%   GATES is a struct array with fields .sites and .gate.
N=2^n;
if strcmpi(flag,'size'),Y=[N N];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_circuit_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=N
    error('um_circuit_apply:Dimension','Input must have 2^n rows.');
end
Y=X;
for layer=1:depth
    if adj,indices=numel(gates):-1:1;else,indices=1:numel(gates);end
    for j=indices
        G=gates(j).gate;
        if adj,G=G';end
        Y=um_gate_apply(Y,G,gates(j).sites,n);
    end
end
end
