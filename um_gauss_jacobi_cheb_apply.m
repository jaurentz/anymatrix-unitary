function Y=um_gauss_jacobi_cheb_apply(X,flag,n,kind)
%UM_GAUSS_JACOBI_CHEB_APPLY FFT-fast exact Gauss-Chebyshev I/II Q and Q'.
%   For ASCENDING Gauss nodes and positive-leading orthonormal
%   polynomials: Q=orthonormal DCT-III*diag((-1)^k) (Chebyshev I);
%                Q=orthonormal DST-I  *diag((-1)^k) (Chebyshev II).
%   The matched DCT/DST kernels use FFT and preserve complex linearity.
%   Both cost O(n log n) per RHS with O(n) transient storage.
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_gauss_jacobi_cheb_apply:Dimension', ...
        'Input must have n rows.');
end
s=(-1).^(0:n-1)';
switch lower(char(flag))
    case 'notransp'
        if kind==1
            Y=um_dct_apply(bsxfun(@times,s,full(X)),'notransp',n,3);
        else
            Y=um_dst_apply(bsxfun(@times,s,full(X)),'notransp',n,1);
        end
    case {'transp','ctransp'}
        if kind==1
            Y=um_dct_apply(full(X),'transp',n,3);
        else
            Y=um_dst_apply(full(X),'transp',n,1);
        end
        Y=bsxfun(@times,s,Y);
    otherwise
        error('um_gauss_jacobi_cheb_apply:Flag', ...
            'Use notransp, transp, ctransp or size.');
end
end
