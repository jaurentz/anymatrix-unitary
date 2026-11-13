function Y=um_scattering_apply(X,flag,inputs,outputs,gates,phase)
%UM_SCATTERING_APPLY Matrix-free disjoint two-port unitary network.
%   Each gate has two input and two output links (ordered, disjoint).
%   INPUTS,OUTPUTS are 2-by-m; GATES is 4-by-m, column-major blocks.
%   The network U=S*diag(PHASE) applies its link phases before scattering.
%   No sparse matrix is assembled. O(n) data and O(n*number_of_rhs)
%   workspace; action cost O(n*number_of_rhs).
n=numel(phase);
if strcmpi(flag,'size'),Y=[n n];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_scattering_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_scattering_apply:Dimension','Input must have n rows.');
end
if size(inputs,1)~=2 || ~isequal(size(inputs),size(outputs)) ...
   || ~isequal(size(gates),[4 size(inputs,2)])
    error('um_scattering_apply:Geometry','Invalid two-port network geometry.');
end
a=inputs(1,:)'; b=inputs(2,:)';
r=outputs(1,:)';s=outputs(2,:)';
g11=gates(1,:)';g21=gates(2,:)';
g12=gates(3,:)';g22=gates(4,:)';
Y=complex(zeros(size(X)));
if adj
    Y(a,:)=conj(g11).*X(r,:)+conj(g21).*X(s,:);
    Y(b,:)=conj(g12).*X(r,:)+conj(g22).*X(s,:);
    Y=conj(phase(:)).*Y;
else
    T=phase(:).*X;
    Y(r,:)=g11.*T(a,:)+g12.*T(b,:);
    Y(s,:)=g21.*T(a,:)+g22.*T(b,:);
end
end
