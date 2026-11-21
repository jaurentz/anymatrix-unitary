function info=um_gauss_moment_diagnostics(nodes,logw,family,alpha)
%UM_GAUSS_MOMENT_DIAGNOSTICS Independent Gaussian quadrature moment audit.
%   INFO compares quadrature (nodes, logarithmic normalized weights)
%   with exact low-degree moments of the NORMALIZED continuous measure.
%   The weights are evaluated independently of the orientation of
%   Golub--Welsch eigenvectors (UM_GAUSS_LOG_WEIGHTS).
%
%   FAMILy='hermite': normalized pi^(-1/2)*exp(-x^2), with
%       m_0=1, m_{2k}=(2k-1)!!/2^k, m_{2k+1}=0.
%   FAMILy='laguerre': normalized Gamma(a)^(-1)*x^(a-1)*exp(-x),
%       a=alpha+1>0; m_k=Gamma(a+k)/Gamma(a)=(a)_k.
%   Gauss N-node quadrature is EXACT for polynomials degree <=2*N-1
%   in mathematical arithmetic. Only orders 0..min(8,2*N-1) checked.
%   These computed residuals diagnose finite-precision accuracy; no
%   claim of a rigorous bound or a fast transform is implied.
%
%   Max-scaled residual divides by max(1,abs(exact moment)).
%   This prevents near-zero odd Hermite moments from spuriously
%   inflating relative errors, while retaining their absolute errors.
%   For extreme alpha causing unrepresentable target/observed moments,
%   diagnostics are explicitly marked as non-finite, not certified.
if nargin<4,alpha=0;end
x=nodes(:);lw=logw(:);
n=numel(x);
if n<1||numel(lw)~=n||any(~isfinite(lw)) ...
        ||any(~isfinite(x))
    error('um_gauss_moment_diagnostics:Input', ...
        'Require finite nodes and finite log weights of equal length.');
end
if ~(ischar(family)||(isstring(family)&&isscalar(family)))
    error('um_gauss_moment_diagnostics:Family','Invalid family.');
end
family=lower(char(family));
if ~any(strcmp(family,{'hermite','laguerre'}))
    error('um_gauss_moment_diagnostics:Family','Unknown Gaussian family.');
end
if strcmp(family,'laguerre')&& ...
        ~(isscalar(alpha)&&isnumeric(alpha)&&isreal(alpha) ...
        &&isfinite(alpha)&&alpha>-1)
    error('um_gauss_moment_diagnostics:Alpha','Alpha must exceed -1.');
end
degree=min(8,2*n-1);
orders=(0:degree)';
w=exp(lw); % tail terms may vanish: record them through log weights
target=zeros(degree+1,1);
measured=zeros(degree+1,1);
target(1)=1;
power=ones(n,1);
for k=0:degree
    if k>0
        power=power.*x;
        if strcmp(family,'laguerre')
            target(k+1)=target(k)*(alpha+k);
        elseif mod(k,2)==0
            target(k+1)=target(k-1)*(k-1)/2;
        end
    end
    measured(k+1)=sum(w.*power);
end
errorAbs=abs(measured-target);
scaled=errorAbs./max(1,abs(target));
allFinite=all(isfinite(target)) && all(isfinite(measured)) ...
    &&all(isfinite(scaled));
if allFinite
    maxScaled=max(scaled);
else
    maxScaled=NaN;
end
info=struct('moment_orders',orders,'exact_moments',target, ...
    'quadrature_moments',measured,'moment_absolute_errors',errorAbs, ...
    'moment_scaled_errors',scaled, ...
    'maximum_scaled_moment_error',maxScaled, ...
    'moment_finite',allFinite, ...
    'moment_degree_exact_in_arithmetic',2*n-1, ...
    'moment_maximum_degree_checked',degree, ...
    'moment_measure',family, ...
    'moment_note',['Independent Gaussian quadrature moment audit ', ...
    '(not transform eigenphase moments and not a proven ', ...
    'floating-point accuracy certificate).']);
end
