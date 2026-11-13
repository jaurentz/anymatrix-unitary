function Y=um_core_apply(X,flag,alpha,beta,steps)
%UM_CORE_APPLY Apply ordered 2-by-2 cores to one or many RHS.
%   'transp' denotes the conjugate transpose, also accepted as 'ctransp'.
n=numel(alpha)+1;
if strcmpi(flag,'size'),Y=[n n];return;end
switch lower(char(flag))
    case 'notransp'
        ordered=steps(end:-1:1);adjoint=false;
    case {'transp','ctransp'}
        ordered=steps;adjoint=true;
    otherwise
        error('um_core_apply:Flag','Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_core_apply:Dimension','Input must have n rows.');
end
Y=X;
for k=ordered
    if k==0
        if adjoint,Y(n,:)=beta*Y(n,:);
        else,Y(n,:)=conj(beta)*Y(n,:);end
    else
        a=alpha(k);rho=sqrt((1-abs(a))*(1+abs(a)));
        if adjoint,G=[a,rho;rho,-conj(a)];
        else,G=[conj(a),rho;rho,-a];end
        Y(k:k+1,:)=G*Y(k:k+1,:);
    end
end
end
