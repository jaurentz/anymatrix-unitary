function Y=um_gauss_jacobi_cheb34_apply(X,flag,n,kind)
%UM_GAUSS_JACOBI_CHEB34_APPLY FFT Gauss-Chebyshev III/IV transform.
%   Exact mathematical Gauss-Jacobi evaluation matrix, nodes ASCENDING:
%       kind=3: (alpha,beta)=(-1/2,+1/2), Chebyshev V_k,
%         theta_j=(2*j-1)*pi/(2*n+1), j=n:-1:1,
%         Q(j,k+1)=2/sqrt(2*n+1)*cos((k+1/2)*theta_j).
%       kind=4: (alpha,beta)=(+1/2,-1/2), Chebyshev W_k,
%         theta_j=2*j*pi/(2*n+1), j=n:-1:1,
%         Q(j,k+1)=2/sqrt(2*n+1)*sin((k+1/2)*theta_j).
%   All columns use positive-leading orthonormal Jacobi polynomials.
%   The roots and positive Gauss weights give Q'*Q=I in exact arithmetic.
%
%   Let M=2*n+1, P=4*M, odd frequencies r=2*k+1, positions
%   s=(2*j-1) [III] or s=2*j [IV]. Then
%       exp(2*pi*i*r*s/P)=exp(i*(k+1/2)*theta_j).
%   Embedding both +/- r Fourier frequencies (or +/- s positions
%   for the transpose) in one length-P IFFT supplies a COMPLEX-LINEAR
%   cosine/sine transform. No dense n-by-n matrix is constructed.
%   FFT length 8*n+4: O(n log n) work, O(n) temporary storage per RHS.
%   Actual floating-point actions are approximate, not exactly isometric.
%
%   FLAG='notransp', 'transp', 'ctransp', or 'size'; Q is real, so
%   transp and ctransp coincide. Complex multi-RHS supported.
%
%   Sources: DLMF 18.3 (Chebyshev III/IV), 3.5 (Gaussian quadrature).
%   https://dlmf.nist.gov/18.3
%   https://dlmf.nist.gov/3.5
if strcmpi(flag,'size')
    Y=[n n]; return
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_gauss_jacobi_cheb34_apply:Dimension', ...
        'Input must be a numeric matrix with n rows.');
end
if ~(kind==3||kind==4)
    error('um_gauss_jacobi_cheb34_apply:Kind', ...
        'Expected third or fourth Chebyshev kind.');
end
flag=lower(char(flag));
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_gauss_jacobi_cheb34_apply:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
m=2*(0:n-1)'+1;
if kind==3
    ell=2*(1:n)'-1;
else
    ell=2*(1:n)';
end
M=2*n+1; P=4*M;
V=complex(zeros(P,size(X,2)));
Z=full(X);
if strcmp(flag,'notransp')
    if kind==3
        V(m+1,:)=Z/2;
        V(P-m+1,:)=Z/2;
    else
        V(m+1,:)=Z/(2i);
        V(P-m+1,:)=-Z/(2i);
    end
    V=ifft(V,[],1);
    Y=flipud(V(ell+1,:))*(2*P/sqrt(M));
else
    Z=flipud(Z);
    if kind==3
        V(ell+1,:)=Z/2;
        V(P-ell+1,:)=Z/2;
    else
        V(ell+1,:)=Z/(2i);
        V(P-ell+1,:)=-Z/(2i);
    end
    V=ifft(V,[],1);
    Y=V(m+1,:)*(2*P/sqrt(M));
end
if isreal(X),Y=real(Y);end
end
