function Y=um_fourier_zernike_apply(X,flag,Q,nr,nphi)
%UM_FOURIER_ZERNIKE_APPLY Fourier synthesis of orthogonal Jacobi radial modes.
%   U=(F_phi^* tensor I_nr)*blockdiag(Q_m), where Q_m are real
%   orthogonal; in MATLAB vectorization, each Q_m acts on column m.
%   Input and output radial slots are the fast-changing indices.
%   The output is NOT samples at the same radius for every m.
%
%   'notransp' : per-mode Q_m, then normalized angular ifft.
%   'ctransp'  : normalized angular fft, then per-mode Q_m'.
%   'transp'   : conjugate of adjoint applied to conjugate input.
%   'size'     : [nr*nphi nr*nphi].
N=nr*nphi;
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_fourier_zernike_apply:Flag','FLAG must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[N N];return;end
if ~any(strcmp(flag,{'notransp','ctransp','transp'}))
    error('um_fourier_zernike_apply:Flag', ...
        'Use notransp, ctransp, transp or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_fourier_zernike_apply:Dimension', ...
        'Input must have nRad*nPhi rows.');
end
if issparse(X),X=full(X);end
transposeOnly=strcmp(flag,'transp');
if transposeOnly,X=conj(X);end
Y=zeros(size(X),'like',complex(0));
for rhs=1:size(X,2)
    C=reshape(X(:,rhs),nr,nphi);
    if strcmp(flag,'notransp')
        for col=1:nphi
            C(:,col)=Q{col}*C(:,col);
        end
        C=sqrt(nphi)*ifft(C,[],2);
    else
        C=fft(C,[],2)/sqrt(nphi);
        for col=1:nphi
            C(:,col)=Q{col}'*C(:,col);
        end
    end
    Y(:,rhs)=C(:);
end
if transposeOnly,Y=conj(Y);end
end
