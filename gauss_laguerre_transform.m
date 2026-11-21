function [U,properties]=gauss_laguerre_transform(n,alpha,output)
%GAUSS_LAGUERRE_TRANSFORM Real orthogonal generalized Gauss-Laguerre Q.
%   Q=GAUSS_LAGUERRE_TRANSFORM(N,ALPHA) returns the positive-leading
%   weighted polynomial evaluation matrix for the normalized measure
%   x^alpha*exp(-x)/Gamma(alpha+1), alpha>-1.
%   Q(j,k)=sqrt(w_j/Gamma(alpha+1))*p_(k-1)(x_j), Q'*Q=I
%   in exact arithmetic. The Jacobi diagonal is 2*k+alpha+1
%   (k=0:N-1), subdiagonal sqrt(k*(k+alpha)) (k=1:N-1).
%
%   'quadrature' returns real nodes, normalized Gauss weights,
%   scaled-log Christoffel weights and orthogonality/recurrence
%   diagnostics. First/last eigenvector endpoint signs are chosen
%   using Sturm oscillation to avoid tail first-component underflow.
%   For L2(dx) Laguerre FUNCTIONS
%      phi_k(x)=x^(alpha/2)*exp(-x/2)*p_k(x)/sqrt(Gamma(alpha+1))
%   we have phi_(k-1)(x_j)=Q(j,k)/d_j, where
%      log(d_j) = .5*(gammaln(alpha+1)+log(w_j/Gamma(alpha+1))
%                      +x_j-alpha*log(x_j)).
%   We expose this as 'log_function_scaling', with no claim that
%   physical function samples are a second orthogonal matrix.
%   This distinction matters near alpha=-1 and in far Gaussian tails.
%
%   'spectral' returns UNKNOWN matrix eigenphases; real Gaussian
%   nodes belong to the Jacobi tridiagonal, not to Q's eigenphases.
%   Full dense O(N^3) Golub-Welsch setup remains; no fast
%   generalized Laguerre operator is claimed.
%   [104] https://doi.org/10.1090/S0025-5718-69-99647-1
%   [106] https://www.cambridge.org/core/journals/acta-numerica/article/abs/fast-algorithms-using-orthogonal-polynomials/4FAD8C7C28EC20EE7583465C1A89AA3D
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(alpha),alpha=0;end
validateattributes(alpha,{'numeric'},{'scalar','real','finite','>',-1});
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('gauss_laguerre_transform:Mode','Invalid output mode.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('gauss_laguerre_transform',n);
    U.alpha=alpha;
    U.notes=['Generalized Gauss-Laguerre nodes are positive Jacobi ', ...
        'eigenvalues, NOT eigenphases of the real orthogonal Q.'];
    return
end
if ~any(strcmp(output,{'matrix','quadrature'}))
    error('gauss_laguerre_transform:Mode', ...
        'Use matrix, quadrature or spectral; no fast operator exists.');
end
k=(0:n-1)';j=(1:n-1)';
d=2*k+alpha+1;e=sqrt(j.*(j+alpha));
if strcmp(output,'matrix')
    U=um_gauss_basis(d,e);
else
    [~,U]=um_gauss_basis(d,e);
    U.family='gauss_laguerre_transform';
    U.alpha=alpha;
    U.log_function_scaling=.5*(gammaln(alpha+1)+ ...
        U.log_normalized_weights+U.nodes-alpha*log(U.nodes));
    U.note=['Physical Laguerre L2(dx) functions are weighted polynomial ', ...
        'samples divided by exp(log_function_scaling); Gaussian ', ...
        'tail weights may require logarithmic representation.'];
end
end
