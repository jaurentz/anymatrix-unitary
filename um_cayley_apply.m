function Y=um_cayley_apply(X,flag,plusFactor,minusFactor,t,n)
%UM_CAYLEY_APPLY Cayley map and adjoint via cached sparse solves.
%   U=(J-itI)*(J+itI)^(-1)=I-2it*(J+itI)^(-1).
if strcmpi(flag,'size'),Y=[n n];return;end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_cayley_apply:Dimension','Input must have n rows.');
end
switch lower(char(flag))
    case 'notransp'
        Y=X-2i*t*(plusFactor\X);
    case {'transp','ctransp'}
        Y=X+2i*t*(minusFactor\X);
    otherwise
        error('um_cayley_apply:Flag','Unknown operator action.');
end
end
