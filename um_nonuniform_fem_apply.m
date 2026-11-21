function Y=um_nonuniform_fem_apply(X,flag,E,L0,alpha,h,Rk,Rm,direct,reverse)
%UM_NONUNIFORM_FEM_APPLY Sparse rational Gauss steps with energy whitening.
%   Physical E*x'=L0*x, normalized U=T*S/T, with
%   T=diag(chol(K),inv(chol(Mass)')). No explicit inverse is formed.
%   Cached generalized LU Cayley factors yield S and S^(-1).
%   S is real, so transpose and conjugate transpose coincide.
n=size(Rk,1);N=2*n;
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_nonuniform_fem_apply:Flag','Invalid operator action.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[N N];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_nonuniform_fem_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_nonuniform_fem_apply:Dimension','Input must have 2N rows.');
end
originalReal=isreal(X);
% To physical canonical x=[q;p]: T^(-1)=[Rk^(-1), Rm'].
Y=[Rk\X(1:n,:);Rm'*X(n+1:end,:)];
if strcmp(flag,'notransp')
    for j=1:numel(alpha)
        Y=direct{j}\((E+alpha(j)*h*L0)*Y);
    end
else
    for j=numel(alpha):-1:1
        Y=reverse{j}\((E-alpha(j)*h*L0)*Y);
    end
end
% Back to energy-normalized coordinates.
Y=[Rk*Y(1:n,:);Rm'\Y(n+1:end,:)];
if originalReal,Y=real(Y);end
end
