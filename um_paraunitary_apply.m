function Y=um_paraunitary_apply(X,flag,L,directions,U0)
%UM_PARAUNITARY_APPLY Circular rank-one delay sections and adjoint.
d=size(directions,1);n=d*L;
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_paraunitary_apply:Flag','Invalid operator flag.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_paraunitary_apply:Dimension','Expected D*L rows.');
end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_paraunitary_apply:Flag','Unknown action flag.');
end
if issparse(X),X=full(X);end
nrhs=size(X,2);
Y=reshape(X,d,L,nrhs);
m=size(directions,2);
if strcmp(flag,'notransp')
    stages=1:m;
    shift=1;
else
    Y=reshape(U0'*reshape(Y,d,[]),d,L,nrhs);
    stages=m:-1:1;
    shift=-1;
end
for t=stages
    v=directions(:,t);
    amplitudes=reshape(v'*reshape(Y,d,[]),L,nrhs);
    difference=circshift(amplitudes,shift,1)-amplitudes;
    correction=reshape(v*reshape(difference,1,[]),d,L,nrhs);
    Y=Y+correction;
end
if strcmp(flag,'notransp')
    Y=reshape(U0*reshape(Y,d,[]),d,L,nrhs);
end
Y=reshape(Y,n,nrhs);
end
