function Y=um_spin_apply(X,flag,n,phase,R,ordering)
%UM_SPIN_APPLY A diagonal Ising phase and commuting one-spin rotations.
%   ORDERING='diag_left' means U=D*R^(tensor n).
%   ORDERING='diag_right' means U=R^(tensor n)*D.
N=2^n;
if strcmpi(flag,'size'),Y=[N N];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_spin_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=N
    error('um_spin_apply:Dimension','Input must have 2^n rows.');
end
Y=X;
if (strcmp(ordering,'diag_left') && adj) || ...
   (strcmp(ordering,'diag_right') && ~adj)
    if adj,Y=conj(phase(:)).*Y;else,Y=phase(:).*Y;end
end
if adj,G=R';else,G=R;end
for j=1:n
    Y=um_gate_apply(Y,G,j,n);
end
if (strcmp(ordering,'diag_left') && ~adj) || ...
   (strcmp(ordering,'diag_right') && adj)
    if adj,Y=conj(phase(:)).*Y;else,Y=phase(:).*Y;end
end
end
