function [U,properties]=periodic_rbf_propagator(n,ell,t,mode)
%PERIODIC_RBF_PROPAGATOR Unitary propagation from a periodic Gaussian kernel.
%   U=PERIODIC_RBF_PROPAGATOR(N,ELL,T) returns exp(1i*T*K), where
%      K(j,k)=exp(-2*sin(pi*(j-k)/N)^2/ELL^2), j,k=0,...,N-1.
%   K is the periodic (exponentiated sine-squared) Gaussian RBF Gram
%   matrix, REAL SYMMETRIC positive definite for finite ELL>0.
%   K is NOT orthogonal/unitary and its polar unitary is the trivial I.
%   The exponential exp(1i*T*K), however, is exactly UNITARY and is
%   generally nontrivial. U is complex symmetric and circulant.
%
%   Positive definiteness follows from
%      exp(-kappa*(1-cos(theta)))
%       =exp(-kappa)*sum_{m=-inf}^inf I_|m|(kappa)*exp(1i*m*theta),
%   where kappa=1/ELL^2 and all Fourier-Bessel coefficients are
%   positive. For equispaced sites, finite eigenvalues are exactly
%      mu_k=sum_{r=0}^{N-1}c_r*exp(-2*pi*1i*r*k/N),
%      c_r=exp(-2*sin(pi*r/N)^2/ELL^2);
%   the eigenphases of U are exp(1i*T*mu_k).
%   No generic RBF Gram eigenvectors or polar projection is asserted.
%
%   MODE='matrix' (default): dense U, N<=2048;
%        'operator': O(N log N) FFT multi-RHS forward/adjoint and
%                    ordinary transpose actions, and 'size' query;
%        'kernel': dense REAL POSITIVE DEFINITE (not unitary) K;
%        'parameters': O(N) Fourier eigenvalues of K and U, kernel
%                      first column, exact finite moment handle;
%        'spectral': finite trace-eigenphase measure, exact in
%                    mathematical arithmetic (unknown limit law).
%   'operator' and 'spectral' allow N up to 1e6.
%   FFT/eigenphase evaluation incurs ordinary floating roundoff.
%
%   This is an application-derived SPECIALIZATION of existing
%   phase_circulant; no new algebraic class or generic fast RBF
%   exponentiation is claimed. For nonuniform nodes the matrix is
%   no longer circulant and this FFT method does not apply.
%
%   [171] C. E. Rasmussen & C. K. I. Williams, Gaussian Processes for
%         Machine Learning, MIT Press (2006), chapter 4, periodic kernels.
%   https://gaussianprocess.org/gpml/chapters/RW4.pdf
properties={'square','unitary','normal','circulant','complex', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'}, ...
    {'scalar','integer','real','positive','finite','<=',1e6});
if nargin<2||isempty(ell),ell=.6;end
if nargin<3||isempty(t),t=.7;end
validateattributes(ell,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(t,{'numeric'}, ...
    {'scalar','real','finite'});
if nargin<4||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('periodic_rbf_propagator:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','operator','kernel','parameters','spectral'}))
    error('periodic_rbf_propagator:Mode', ...
        'Use matrix, operator, kernel, parameters or spectral.');
end
if any(strcmp(mode,{'matrix','kernel'})) && n>2048
    error('periodic_rbf_propagator:DenseOrder', ...
        'Dense matrix/kernel modes need N<=2048.');
end
r=(0:n-1)';
c=exp(-2*(sin(pi*r/n)/ell).^2);
mu=real(fft(c));
% Exact K is symmetric real circulant, giving mu(k)=mu(-k).
% Enforce conjugate-pair symmetry within roundoff (not clipping
% finite Gaussian kernel eigenvalues or applying a polar factor).
mirror=[1;(n:-1:2)'];
mu=.5*(mu+mu(mirror));
z=exp(1i*t*mu);
switch mode
    case 'operator'
        U=@(X,flag)um_periodic_rbf_apply(X,flag,z);
    case {'matrix','kernel'}
        if strcmp(mode,'matrix'),column=ifft(z);
        else,column=c;end
        indices=mod((0:n-1)'-(0:n-1),n)+1;
        U=column(indices);
    case 'parameters'
        U=struct('kernel_first_column',c, ...
            'kernel_eigenvalues',mu, ...
            'propagator_eigenvalues',z, ...
            'moments',@(k)rbf_finite_moments(z,k), ...
            'dimension',n,'lengthscale',ell,'time',t, ...
            'kernel','periodic exponentiated sine squared', ...
            'reference','Rasmussen--Williams (2006), ch. 4', ...
            'notes',['The real kernel K is positive definite, ', ...
                'NOT unitary; U=exp(1i*t*K) is unitary. ', ...
                'Kernel eigenvalues mu are real and are NOT ', ...
                'the eigenphases of U; the finite phase law ', ...
                'is exp(1i*t*mu). Large kernel spectra can ', ...
                'be near machine precision.']);
    case 'spectral'
        U=unitary_spectral_reference('periodic_rbf_propagator',n);
        U.finite.kind='exact';
        U.finite.eigenvalues=@()z;
        U.finite.moments=@(k)rbf_finite_moments(z,k);
        U.reference='Rasmussen--Williams (2006), chapter 4';
        U.notes=['Exactly known finite circulant RBF propagator ', ...
            'trace phases exp(1i*t*mu); kernel K is SPD but ', ...
            'NOT unitary. No universal limiting law claimed.'];
        U.lengthscale=ell;U.time=t;
end
end

function m=rbf_finite_moments(z,k)
validateattributes(k,{'numeric'}, ...
    {'real','integer','finite'});
m=reshape(arrayfun(@(q)mean(z.^q),k(:)),size(k));
end
