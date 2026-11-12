function [U,order] = js_core_product(alpha,beta,representation,twists)
%JS_CORE_PRODUCT Finite unitary multiplication matrices from common cores.
%   U = JS_CORE_PRODUCT(ALPHA,BETA,REPRESENTATION,TWISTS) uses
%   alpha_0,...,alpha_(n-2), |alpha_j|<1, and terminal beta in T.
%   Theta(a) = [conj(a),rho;rho,-a], rho=sqrt(1-|a|^2).
%
%   REPRESENTATION = 'hessenberg', 'inverse_hessenberg', 'cmv', 'twisted'.
%   In the twisted case TWISTS is a logical vector of length n-2:
%   starting from core index 1, for each new core k=2,...,n-1,
%   1 prepends core k, 0 appends it. A default alternating pattern is
%   supplied. All factors occur once. The terminal diagonal D_n(beta)
%   is placed at the right for the generic snake form.
%
%   CMV uses its exact parity-dependent L*M factorization, while inverse
%   Hessenberg is D*G_(n-1)*...*G_1. All variants share finite Szego
%   quadrature eigenvalues and e1 weights (with the same beta).
if nargin < 3 || isempty(representation), representation='cmv'; end
if nargin < 4, twists=[]; end
if ~isempty(alpha)
    validateattributes(alpha,{'numeric'},{'vector','finite'});
elseif ~isnumeric(alpha)
    error('js_core_product:Alpha','ALPHA must be numeric.');
end
alpha=alpha(:);
n=numel(alpha)+1;
if any(abs(alpha)>=1)
    error('js_core_product:Alpha','Interior Verblunsky coefficients must have modulus < 1.');
end
validateattributes(beta,{'numeric'},{'scalar','finite'});
if abs(abs(beta)-1) > 128*eps
    error('js_core_product:Beta','Closure beta must be unimodular.');
end
beta=beta/abs(beta);
representation=lower(char(representation));
order=[];
switch representation
    case 'hessenberg'
        order=1:n-1;
        U=speye(n);
        for k=order
            U=apply_right(U,k,alpha(k));
        end
        U(:,n)=conj(beta)*U(:,n);
    case {'inverse_hessenberg','inverse'}
        order=n-1:-1:1;
        U=speye(n);
        U(n,n)=conj(beta);
        for k=order
            U=apply_right(U,k,alpha(k));
        end
    case 'cmv'
        L=speye(n); M=speye(n);
        for k=1:2:n-1
            L=apply_right(L,k,alpha(k));
        end
        for k=2:2:n-1
            M=apply_right(M,k,alpha(k));
        end
        if mod(n,2)==1
            L(n,n)=conj(beta);
        else
            M(n,n)=conj(beta);
        end
        U=L*M;
        order=[1:2:n-1, 2:2:n-1];
    case 'twisted'
        if isempty(twists)
            twists=mod((1:max(n-2,0))',2)==1;
        end
        validateattributes(twists,{'numeric','logical'},{'vector'});
        if numel(twists)~=max(0,n-2) || any(twists(:)~=0 & twists(:)~=1)
            error('js_core_product:Twists', ...
                'TWISTS must have n-2 binary entries.');
        end
        order=1:min(n-1,1);
        for k=2:n-1
            if twists(k-1)
                order=[k,order]; %#ok<AGROW>
            else
                order=[order,k]; %#ok<AGROW>
            end
        end
        U=speye(n);
        for k=order
            U=apply_right(U,k,alpha(k));
        end
        U(:,n)=conj(beta)*U(:,n);
    otherwise
        error('js_core_product:Representation','Unknown representation %s.',representation);
end
end

function U=apply_right(U,k,a)
% Product by an embedded elementary core without forming an n-by-n core.
rho=sqrt((1-abs(a))*(1+abs(a)));
U(:,k:k+1)=U(:,k:k+1)*[conj(a),rho;rho,-a];
end
