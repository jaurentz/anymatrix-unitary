function Y=um_gyroscopic_oscillator_apply(X,flag,P,Q,T)
%UM_GYROSCOPIC_OSCILLATOR_APPLY Sparse rational normalized actions.
% U=T*(Q\P)/T. Its adjoint is the inverse T*(P\Q)/T.
n=size(P,1);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_gyroscopic_oscillator_apply:Flag','Invalid flag.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_gyroscopic_oscillator_apply:Flag','Unknown action.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_gyroscopic_oscillator_apply:Dimension','Expected 2N rows.');
end
if strcmp(flag,'notransp')
    Y=T*(Q\(P*(T\X)));
else
    Y=T*(P\(Q*(T\X)));
end
end
