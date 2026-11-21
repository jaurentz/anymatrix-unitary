function Y=um_tensor_fem_apply(X,flag,grid,z)
%UM_TENSOR_FEM_APPLY Fast separable FE normalized orthogonal propagator.
%   The real step is [C,S;-S,C] with modal C=Re(z), S=-Im(z).
%   The d-dimensional DST is its own inverse. For complex RHS,
%   'transp' and 'ctransp' agree because the unitary is real.
N=prod(grid);
if strcmpi(flag,'size'),Y=[2*N,2*N];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_tensor_fem_apply:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=2*N
    error('um_tensor_fem_apply:Dimension', ...
        'Input must have 2*prod(GRID) rows.');
end
q=um_tensor_dst_apply(X(1:N,:),grid);
p=um_tensor_dst_apply(X(N+1:end,:),grid);
a=real(z(:));b=-imag(z(:));
if adj,b=-b;end
Y=[um_tensor_dst_apply(a.*q+b.*p,grid); ...
   um_tensor_dst_apply(-b.*q+a.*p,grid)];
end
