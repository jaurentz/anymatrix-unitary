function Y=um_kick_tridiag_apply(X,flag,phase,d,u,t)
%UM_KICK_TRIDIAG_APPLY Diagonal phase times tridiagonal exponential.
%   Operator U=diag(phase)*exp(-1i*t*H), with Hermitian tridiagonal H
%   represented by diagonal d and upper diagonal u.
%   'transp' and 'ctransp' are both the complex conjugate adjoint.
n=numel(phase);
if strcmpi(flag,'size'),Y=[n n];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_kick_tridiag_apply:Flag', ...
          'Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_kick_tridiag_apply:Dimension','Input must have n rows.');
end
if adj
    Y=um_tridiag_exp_apply(conj(phase(:)).*X,d,u,-t);
else
    Y=phase(:).*um_tridiag_exp_apply(X,d,u,t);
end
end
