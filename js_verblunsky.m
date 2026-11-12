function alpha = js_verblunsky(mu,count)
%JS_VERBLUNSKY Return interior OPUC coefficients for a measure descriptor.
%   This dispatcher is the extension point for new Szego measures.
%   The coefficient sequence here does NOT include a finite unitary closure.
assert(isstruct(mu) && isfield(mu,'family'),'Expected a measure descriptor.');
switch mu.family
    case 'jacobi_szego'
        alpha=js_verblunsky_jacobi(count,mu.a,mu.b);
    case 'bernstein_szego'
        validateattributes(count,{'numeric'},{'scalar','integer','nonnegative','finite'});
        alpha=zeros(count,1);
        if count>0,alpha(1)=mu.alpha0;end
    otherwise
        error('js_verblunsky:Family','No coefficient generator for %s.',mu.family);
end
end
