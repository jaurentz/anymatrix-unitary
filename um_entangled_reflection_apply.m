function Y=um_entangled_reflection_apply(X,flag,v,lambda,n)
%UM_ENTANGLED_REFLECTION_APPLY Rank-one phase reflection action.
%   Spectral projector P=v*v' gives U=I+(lambda-1)*P.
if strcmpi(flag,'size'),Y=[n n];return;end
if any(strcmpi(flag,{'transp','ctransp'}))
    c=conj(lambda)-1;
elseif strcmpi(flag,'notransp')
    c=lambda-1;
else
    error('um_entangled_reflection_apply:Flag','Invalid operator action.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_entangled_reflection_apply:Dimension','Input must have d^2 rows.');
end
Y=X+c*v*(v'*X);
end
