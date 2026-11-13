function alpha = js_verblunsky(mu,count)
%JS_VERBLUNSKY Interior Verblunsky coefficients for probability measures.
%   ALPHA contains coefficients alpha_0,...,alpha_(COUNT-1), without
%   the finite terminal beta. See SZEGO_MATRIX for finite unitary closure.
assert(isstruct(mu) && isfield(mu,'family'),'Expected a measure descriptor.');
validateattributes(count,{'numeric'},{'scalar','integer','nonnegative','finite'});
switch mu.family
    case 'jacobi_szego'
        alpha=js_verblunsky_jacobi(count,mu.a,mu.b);
    case 'bernstein_szego'
        alpha=zeros(count,1);
        if count>0,alpha(1)=mu.alpha0;end
    case 'bernstein_finite'
        alpha=zeros(count,1);
        m=min(count,numel(mu.coefficients));
        alpha(1:m)=mu.coefficients(1:m);
    case 'rogers_szego'
        k=(0:count-1)';
        alpha=(-1).^k.*mu.q.^((k+1)/2);
    case 'circular_jacobi'
        alpha=zeros(count,1);
        b=mu.a+1i*mu.eta;
        % Pochhammer ratio: -(b)_(k+1)/(1+conj(b))_(k+1).
        value=1;
        for k=0:count-1
            value=value*(b+k)/(1+conj(b)+k);
            alpha(k+1)=-value;
        end
    case 'geronimus'
        alpha=repmat(mu.alpha,count,1);
    case 'periodic'
        alpha=mu.pattern(mod((0:count-1)',mu.period)+1);
    otherwise
        error('js_verblunsky:Family','No coefficient generator for %s.',mu.family);
end
if any(~isfinite(alpha)) || any(abs(alpha)>=1)
    error('js_verblunsky:Validity', ...
        'The requested coefficients are not strictly within the unit disk.');
end
end
