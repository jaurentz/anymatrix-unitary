function [U,properties]=gauss_hermite_transform(n,output)
%GAUSS_HERMITE_TRANSFORM Real orthogonal Gauss-Hermite evaluation.
%   Q=GAUSS_HERMITE_TRANSFORM(N) maps positive-leading orthonormal
%   physicists' Hermite POLYNOMIAL coefficients to quadrature-weighted
%   values, Q(j,k)=sqrt(w_j/sqrt(pi))*p_(k-1)(x_j). N-point Gaussian
%   exactness implies Q'*Q=I mathematically.
%   Jacobi diagonal 0, offdiagonal sqrt(k/2), k=1,...,N-1.
%
%   'quadrature' also audits the exact normalized Hermite measure
%   moments m_(2k)=(2k-1)!!/2^k, m_(2k+1)=0, for degrees up to
%   min(8,2*n-1), independently of eigenvector sign/Gram residuals.
%   Reported moment errors are measured, not rigorous upper bounds.
%   'quadrature' returns a diagnostic struct with sorted real
%   Hermite nodes, normalized weights and LOG normalized weights
%   from rescaled Christoffel polynomials. Tiny tail weights can
%   underflow to zero as doubles while their logarithms remain finite.
%   Oriented Jacobi eigenvectors use a Sturm first/last endpoint
%   sign relation, avoiding loss of sign when first components
%   underflow. Gram and Jacobi diagonalization residuals are reported.
%
%   IMPORTANT FUNCTION versus POLYNOMIAL normalization: let
%      phi_k(x)=pi^(-1/4)*exp(-x^2/2)*p_k(x)
%   be the physically L2(dx)-orthonormal Hermite FUNCTION. Then
%      phi_(k-1)(x_j)=Q(j,k)/d_j,
%      log(d_j)=.25*log(pi)+.5*log(w_j/sqrt(pi))+.5*x_j^2.
%   The positive diagonal scaling is returned as
%   'log_function_scaling' without exponentiating large factors.
%   It does not define an additional unitary matrix.
%
%   'spectral' returns UNKNOWN matrix eigenphases. Gaussian nodes
%   are the real Jacobi eigenvalues, NOT eigenphases of Q.
%   Full Q still requires dense O(N^3) eigensolver setup and
%   O(N^2) storage/application. No fast operator is claimed.
%   [104] https://doi.org/10.1090/S0025-5718-69-99647-1
%   [109] https://arxiv.org/abs/2604.02041
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('gauss_hermite_transform:Mode','Invalid output mode.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('gauss_hermite_transform',n);
    U.notes=['Gauss-Hermite quadrature nodes are Jacobi eigenvalues, ' ...
        'NOT eigenphases of the orthogonal transform.'];
    return
end
if ~any(strcmp(output,{'matrix','quadrature'}))
    error('gauss_hermite_transform:Mode', ...
        'Use matrix, quadrature or spectral; no fast operator exists.');
end
d=zeros(n,1);e=sqrt((1:n-1)'/2);
if strcmp(output,'matrix')
    U=um_gauss_basis(d,e);
else
    [~,U]=um_gauss_basis(d,e);
    U.family='gauss_hermite_transform';
    U.moment_diagnostics=um_gauss_moment_diagnostics( ...
        U.nodes,U.log_normalized_weights,'hermite');
    U.log_function_scaling=.25*log(pi)+ ...
        .5*U.log_normalized_weights+.5*U.nodes.^2;
    U.note=['The quadrature operator Q is orthogonal; physical Hermite ', ...
        'function samples are Q divided rowwise by exp(log_function_scaling). ', ...
        'Exponential tails may underflow in weights but remain in logarithmic form.'];
end
end
