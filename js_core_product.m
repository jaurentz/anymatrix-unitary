function [U,order]=js_core_product(alpha,beta,representation,twists,output)
%JS_CORE_PRODUCT Shared explicit and implicit finite Verblunsky products.
%   OUTPUT='operator' gives U(X,'notransp'/'transp'/'size').
if nargin<3 || isempty(representation),representation='cmv';end
if nargin<4,twists=[];end
if nargin<5,output='matrix';end
if ~isempty(alpha)
    validateattributes(alpha,{'numeric'},{'vector','finite'});
elseif ~isnumeric(alpha)
    error('js_core_product:Alpha','ALPHA must be numeric.');
end
alpha=alpha(:);n=numel(alpha)+1;
if any(abs(alpha)>=1)
    error('js_core_product:Alpha','Interior Verblunsky coefficients must have modulus < 1.');
end
validateattributes(beta,{'numeric'},{'scalar','finite'});
if abs(abs(beta)-1)>128*eps
    error('js_core_product:Beta','Closure beta must be unimodular.');
end
beta=beta/abs(beta);
[order,steps]=um_core_order(n,representation,twists);
if um_is_operator(output)
    U=@(X,flag)um_core_apply(X,flag,alpha,beta,steps);
    return
end
U=speye(n);
for k=steps
    if k==0
        U(:,n)=conj(beta)*U(:,n);
    else
        a=alpha(k);rho=sqrt((1-abs(a))*(1+abs(a)));
        U(:,k:k+1)=U(:,k:k+1)*[conj(a),rho;rho,-a];
    end
end
end
