function Y=um_hessenberg_pencil_apply(X,flag,alphaA,betaA,phaseA,alphaB,betaB,steps)
%UM_HESSENBERG_PENCIL_APPLY Exact O(n) core-factor quotient actions.
%   A=diag(phaseA)*H(alphaA,betaA), B=H(alphaB,betaB), U=B'*A.
%   Inputs are unitary core parameter arrays (including reducible |a|=1).
%   'transp' and 'ctransp' both mean conjugate transpose.
n=numel(alphaA)+1;
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_hessenberg_pencil_apply:Dimension', ...
        'Input must have n rows.');
end
switch lower(char(flag))
    case 'notransp'
        Y=um_core_apply(X,'notransp',alphaA,betaA,steps);
        Y=bsxfun(@times,phaseA,Y);
        Y=um_core_apply(Y,'ctransp',alphaB,betaB,steps);
    case {'transp','ctransp'}
        Y=um_core_apply(X,'notransp',alphaB,betaB,steps);
        Y=bsxfun(@times,conj(phaseA),Y);
        Y=um_core_apply(Y,'ctransp',alphaA,betaA,steps);
    otherwise
        error('um_hessenberg_pencil_apply:Flag', ...
            'Use notransp, transp, ctransp or size.');
end
end
