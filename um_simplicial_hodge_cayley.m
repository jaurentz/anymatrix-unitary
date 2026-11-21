function Y=um_simplicial_hodge_cayley(X,flag,N,steps,Ip,Im,plusFactor,minusFactor)
%UM_SIMPLICIAL_HODGE_CAYLEY Sparse unitary rational Hodge-Dirac action.
%   The single Crank–Nicolson step C=(I+i*h*H/2)^(-1)
%                                     *(I-i*h*H/2)
%   is exactly unitary in mathematical arithmetic for H real
%   symmetric, and C^T=C. However C^steps is NOT exp(-itH):
%   for fixed H and t the approximation error is O(steps^-2).
%   Sparse LU factors for Ip/Im are CACHED by the caller.
%   Ordinary transpose C^T=C; adjoint C'=C^(-1).
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_simplicial_hodge_cayley:Flag','FLAG must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[N N];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_simplicial_hodge_cayley:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_simplicial_hodge_cayley:Dimension', ...
        'Input must have cochain dimension rows.');
end
if issparse(X),X=full(X);end
Y=X;
if strcmp(flag,'ctransp')
    for k=1:steps,Y=minusFactor\(Ip*Y);end
else
    % The complex symmetric C has the SAME ordinary transpose.
    for k=1:steps,Y=plusFactor\(Im*Y);end
end
end
