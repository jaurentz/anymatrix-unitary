function Y=um_reflection_apply(X,flag,n,marked)
%UM_REFLECTION_APPLY Grover oracle and rank-one diffusion; exact adjoint.
if strcmpi(flag,'size'),Y=[n n];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_reflection_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_reflection_apply:Dimension','Input must have n rows.');
end
Y=X;
if ~adj,Y(1:marked,:)=-Y(1:marked,:);end
Y=2*mean(Y,1)-Y;
if adj,Y(1:marked,:)=-Y(1:marked,:);end
end
