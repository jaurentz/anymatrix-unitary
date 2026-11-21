function [U,properties]=gauss_jacobi_transform(n,alpha,beta,output)
%GAUSS_JACOBI_TRANSFORM Orthogonal Gauss-Jacobi weighted evaluation.
%   U=GAUSS_JACOBI_TRANSFORM(N,ALPHA,BETA) is the N-by-N orthogonal
%   matrix sqrt(w_j/mu0)*p_{k-1}(x_j), where p_k are orthonormal
%   with respect to the normalized Jacobi measure proportional to
%   (1-x)^alpha*(1+x)^beta on [-1,1] and the (x_j,w_j) are
%   N-point Gauss-Jacobi data. ALPHA,BETA > -1 (default 0,0).
%   ALPHA=BETA=0: Legendre; =-1/2: Chebyshev I; =1/2: Chebyshev II.
%   ALPHA=BETA=lambda-1/2: Gegenbauer family.
%   Quadrature nodes are the Jacobi tridiagonal eigenvalues, NOT
%   the unit-circle eigenvalues of U.
%   OUTPUT='matrix' (default) preserves the exact finite Golub-Welsch
%   matrix for ALL alpha,beta>-1. OUTPUT='operator' supplies a GENUINE
%   O(N log N), O(N) workspace FFT action for the two classical,
%   exactly reducible Chebyshev specializations ONLY:
%      alpha=beta=-1/2: Q=DCT-III * diag((-1)^k) (first kind);
%      alpha=beta= 1/2: Q=DST-I   * diag((-1)^k) (second kind).
%   These signs follow ASCENDING Gauss node ordering and POSITIVE-
%   leading orthonormal Jacobi polynomials. Both forward/adjoint actions
%   are strictly the same mathematical orthogonal matrix as default.
%   General alpha,beta (including associated-Legendre m,m) explicitly
%   REJECT 'operator' pending a verified fast connection algorithm.
%   OUTPUT='spectral' reports UNKNOWN transform eigenphases: the real
%   quadrature nodes are NOT the spectrum of Q.
%   Generic dense Golub-Welsch setup remains O(N^3).
%   [104] https://doi.org/10.1090/S0025-5718-69-99647-1
%   [105] https://doi.org/10.1093/imanum/drv060
%   [106] https://www.cambridge.org/core/journals/acta-numerica/article/abs/fast-algorithms-using-orthogonal-polynomials/4FAD8C7C28EC20EE7583465C1A89AA3D
%   [107] https://arxiv.org/abs/1705.05448
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(alpha),alpha=0;end
if nargin<3||isempty(beta),beta=0;end
if nargin<4||isempty(output),output='matrix';end
validateattributes(alpha,{'numeric'},{'scalar','real','finite','>',-1});
validateattributes(beta,{'numeric'},{'scalar','real','finite','>',-1});
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('gauss_jacobi_transform:Mode','OUTPUT must be scalar text.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('gauss_jacobi_transform',n);
    U.reference='Golub--Welsch Gaussian quadrature and fast Chebyshev FFT';
    U.notes=['The Gauss-Jacobi nodes are Jacobi-matrix real eigenvalues, ', ...
        'NOT the generally unknown eigenphases of Q.'];
    return
end
if any(strcmp(output,{'operator','implicit'}))
    if alpha==-.5 && beta==-.5
        kind=1;
    elseif alpha==.5 && beta==.5
        kind=2;
    else
        error('gauss_jacobi_transform:FastSubset', ...
            ['Fast FFT mode currently covers only alpha=beta=+/-1/2; ', ...
             'no hidden dense general Jacobi operator is supplied.']);
    end
    U=@(X,flag)um_gauss_jacobi_cheb_apply(X,flag,n,kind);
    return
end
if ~any(strcmp(output,{'matrix','explicit'}))
    error('gauss_jacobi_transform:Mode','Unsupported output mode.');
end
k=(0:n-1)';
d=zeros(n,1);
d(1)=(beta-alpha)/(alpha+beta+2);
if n>1
    t=k(2:end);
    d(2:end)=(beta^2-alpha^2)./ ...
        ((2*t+alpha+beta).*(2*t+alpha+beta+2));
end
e=zeros(max(n-1,0),1);
if n>1
    % Special k=1 variance formula avoids a removable 0/0 at
    % alpha+beta=-1 and is well behaved near that line.
    e(1)=2/(alpha+beta+2)*sqrt((1+alpha)*(1+beta)/(alpha+beta+3));
end
if n>2
    t=(2:n-1)';
    e(2:end)=2./(2*t+alpha+beta).* ...
        sqrt(t.*(t+alpha).*(t+beta).*(t+alpha+beta)./ ...
            ((2*t+alpha+beta-1).*(2*t+alpha+beta+1)));
end
U=um_gauss_basis(d,e);
end
