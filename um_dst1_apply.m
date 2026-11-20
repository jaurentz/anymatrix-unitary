function Y=um_dst1_apply(X)
%UM_DST1_APPLY Orthonormal discrete sine transform of type I along rows.
%   Uses a length 2*(n+1) FFT odd extension, with no toolboxes.
%   Self-inverse, including for complex arrays with multiple columns.
if ~isnumeric(X) || ndims(X)>2 || size(X,1)<1
    error('um_dst1_apply:Dimension','Input must have at least one row.');
end
if issparse(X), X=full(X); end
n=size(X,1);
Z=[zeros(1,size(X,2));X;zeros(1,size(X,2));-flipud(X)];
F=fft(Z,[],1);
Y=(1i/2)*sqrt(2/(n+1))*F(2:n+1,:);
end
