function Y=um_geometric_oscillator_apply(X,flag,z)
%UM_GEOMETRIC_OSCILLATOR_APPLY O(n log n) realified unitary oscillator step.
%   Z stores positive-frequency eigenvalues in Fourier mode order.
%   U=[C S;-S C], C=ifft*diag(real(Z))*fft, S=ifft*diag(-imag(Z))*fft.
%   This action supports complex input and multiple right-hand sides.
n=numel(z);
if strcmpi(flag,'size'), Y=[2*n,2*n]; return; end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_geometric_oscillator_apply:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=2*n
    error('um_geometric_oscillator_apply:Dimension', ...
        'Input must have 2*n rows.');
end
if issparse(X), X=full(X); end
a=real(z(:)); b=-imag(z(:));
if adj, b=-b; end
q=fft(X(1:n,:),[],1);
p=fft(X(n+1:end,:),[],1);
Y=[ifft(a.*q+b.*p,[],1);ifft(-b.*q+a.*p,[],1)];
end
