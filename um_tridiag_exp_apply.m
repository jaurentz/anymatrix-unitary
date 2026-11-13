function Y=um_tridiag_exp_apply(X,d,u,t)
%UM_TRIDIAG_EXP_APPLY Matrix-free exponential of Hermitian tridiagonal H.
%   Y=UM_TRIDIAG_EXP_APPLY(X,D,U,T) returns exp(-1i*T*H)*X,
%   where H has real diagonal D, upper diagonal U and lower
%   diagonal conj(U). No dense eigendecomposition or EXPM is formed.
%
%   Uses the Chebyshev-Bessel series on a Gershgorin spectral interval,
%   with a numerically selected tail cutoff. The kernel stores O(n)
%   coefficients and O(n*number_of_rhs) working arrays. Computation
%   scales with the spectral radius times |T|; it is not uniformly
%   O(n) arithmetic for an arbitrary time or spectral range.
%
%   The cutoff estimate checks a generously extended Bessel sequence;
%   numerical truncation error is typically near roundoff in tested
%   regimes. It is NOT a certified arbitrary-precision expmv routine.
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
