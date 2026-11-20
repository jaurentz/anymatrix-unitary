function Y=um_geometric_wave_apply(X,flag,z)
%UM_GEOMETRIC_WAVE_APPLY DST-I matrix-free action for finite-element wave steps.
%   Z contains the negative-frequency phases ordered by Dirichlet sine modes.
n=numel(z);
if strcmpi(flag,'size'),Y=[2*n,2*n];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_geometric_wave_apply:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=2*n
    error('um_geometric_wave_apply:Dimension','Input must have 2*n rows.');
end
if issparse(X), X=full(X); end
q=um_dst1_apply(X(1:n,:));
p=um_dst1_apply(X(n+1:end,:));
a=real(z(:));b=-imag(z(:));
if adj,b=-b;end
Y=[um_dst1_apply(a.*q+b.*p);um_dst1_apply(-b.*q+a.*p)];
end
