function Y=um_periodic_rbf_apply(X,flag,z)
%UM_PERIODIC_RBF_APPLY FFT action for symmetric complex RBF propagator.
%   U=F'*diag(z)*F with z(-k)=z(k) because the kernel is real
%   symmetric circulant. Hence U.'=U, whereas U' generally differs.
n=numel(z);
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_periodic_rbf_apply:Flag','FLAG must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_periodic_rbf_apply:Flag', ...
        'Use notransp, transp, ctransp, or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_periodic_rbf_apply:Dimension', ...
        'Input must have N rows.');
end
if strcmp(flag,'ctransp'),p=conj(z(:));
else,p=z(:);end
Y=um_fft_apply(X,'notransp','circulant',p,[]);
end
