function Y=um_tridiag_exp_apply(X,d,u,t,method)
%UM_TRIDIAG_EXP_APPLY Hermitian-tridiagonal exponential times columns.
%   Y=UM_TRIDIAG_EXP_APPLY(X,D,U,T) computes exp(-1i*T*H)*X for
%   sparse Hermitian tridiagonal H with diagonal D and superdiagonal U.
%
%   On MATLAB R2023b+ the default calls MathWorks EXPMV, one column
%   at a time, with a sparse H. Unlike EXPM it does not form exp(H).
%   On older MATLAB (repository minimum R2020b), a Chebyshev-Bessel
%   exponential action provides a compatible fallback.
%   Internal optional METHOD='chebyshev' tests that fallback explicitly,
%   and METHOD='expmv' requests the MathWorks implementation.
%
%   Neither action is an exactly unitary finite matrix operation in
%   floating point; the mathematical propagator is unitary.
d=d(:);
u=u(:);
n=numel(d);
if n<1 || numel(u)~=n-1
    error('um_tridiag_exp_apply:Shape','Need n real diagonal and n-1 upper entries.');
end
if ~isreal(d) || any(~isfinite(d)) || any(~isfinite(u)) ...
        || ~isscalar(t) || ~isreal(t) || ~isfinite(t)
    error('um_tridiag_exp_apply:Parameters','Invalid Hermitian tridiagonal or time.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_tridiag_exp_apply:Dimension','Input must have n rows.');
end
if nargin<5 || isempty(method),method='auto';end
if ~(ischar(method) || (isstring(method) && isscalar(method)))
    error('um_tridiag_exp_apply:Method','Method must be auto, expmv or chebyshev.');
end
method=lower(char(method));
if ~any(strcmp(method,{'auto','expmv','chebyshev'}))
    error('um_tridiag_exp_apply:Method','Unknown exponential-action method.');
end
hasExpmv=exist('expmv','file')~=0 || exist('expmv','builtin')~=0;
useExpmv=strcmp(method,'expmv') || (strcmp(method,'auto') && hasExpmv);
if useExpmv
    if ~hasExpmv
        error('um_tridiag_exp_apply:Availability', ...
              'EXPMV requires MATLAB R2023b or later.');
    end
    if t==0,Y=X;return;end
    H=um_sparse_hermitian_tridiag(d,u);
    A=-1i*H;
    Y=complex(zeros(size(X)));
    for j=1:size(X,2)
        Y(:,j)=expmv(A,X(:,j),t);
    end
    return
end
% Compatible R2020b fallback: Chebyshev-Bessel on a spectral enclosure.
c=(max(d)+min(d))/2;
a=max(abs(d-c));
if ~isempty(u),a=a+2*max(abs(u));end
if a==0 || t==0
    Y=exp(-1i*t*c)*X;
    return
end
x=abs(t)*a;
if ~isfinite(x) || x>5e4
    error('um_tridiag_exp_apply:Scale', ...
          'Spectral time range too large for the Chebyshev-Bessel action.');
end
% The enlarged range is beyond the Bessel turning point k ~= x.
m=ceil(x+14*sqrt(x+1)+60);
b=besselj((0:m)',t*a);
tail=[flipud(cumsum(flipud(abs(b))));0];
kmax=find(2*tail(2:end)<5e-14,1,'first')-1;
if isempty(kmax),kmax=m;end
Y=b(1)*X;
if kmax>0
    t0=X;
    t1=scaled_action(X,d,u,c,a);
    for k=1:kmax
        if k>1
            t2=2*scaled_action(t1,d,u,c,a)-t0;
            t0=t1;
            t1=t2;
        end
        Y=Y+2*((-1i)^k)*b(k+1)*t1;
    end
end
Y=exp(-1i*t*c)*Y;
end

function Y=scaled_action(X,d,u,c,a)
Y=((d-c)/a).*X;
if ~isempty(u)
    Y(1:end-1,:)=Y(1:end-1,:)+(u/a).*X(2:end,:);
    Y(2:end,:)=Y(2:end,:)+(conj(u)/a).*X(1:end-1,:);
end
end
