function [U,properties]=toda_cayley(n,tau,gamma,a,b,output)
%TODA_CAYLEY Isospectral finite Jacobi Toda flow and its unitary Cayley image.
%   U=TODA_CAYLEY(N,TAU,GAMMA,A,B,OUTPUT) starts from an N-by-N real
%   symmetric irreducible Jacobi matrix J0. By default J0 is the
%   Gauss-Jacobi recurrence for weight (1-x)^A*(1+x)^B, A,B>-1,
%   scaled to [-2,2]; A=B=0 by default. Alternatively A may be a
%   struct with fields 'diagonal' (N values) and 'subdiagonal' (N-1
%   strictly positive values), with B empty.
%
%       exp(TAU*J0) = Q(TAU)*R(TAU), diag(R)>0,
%       J(TAU) = Q(TAU)'*J0*Q(TAU),
%       U = (J(TAU)-1i*GAMMA*I)/(J(TAU)+1i*GAMMA*I).
%
%   The exact continuous finite nonperiodic Toda QR solution is a
%   symmetric tridiagonal J(TAU) with positive subdiagonal. U is
%   complex symmetric and unitary. Its fixed eigenvalues are
%   (x_j-1i*GAMMA)/(x_j+1i*GAMMA), x_j=eig(J0).
%   The e1 CYCLIC weights vary as w_j(TAU) proportional to
%   w_j(0)*exp(2*TAU*x_j), unlike the fixed equal-weight trace measure.
%
%   OUTPUT='matrix' (default), 'operator', 'spectral', 'jacobi', or
%   'parameters'. The operator caches sparse shifted LU solves, but
%   generating J(TAU) by a dense QR factorization has O(N^3) setup.
%   'jacobi' returns a real symmetric tridiagonal, not a unitary matrix.
%   The shifted exponential avoids overflow but the QR trajectory can
%   lose accuracy when exponential weights become ill-conditioned.
%
%   Symes, Physica D 4 (1982), 275--280:
%   https://doi.org/10.1016/0167-2789(82)90069-0
%   Aurentz et al., ETNA 46 (2017), 447--459:
%   https://www.maths.tcd.ie/EMIS/journals/ETNA/volumes/2011-2020/vol46/abstract_vol46_pp447-459.html?pages=447-459&vol=46
properties={'square','unitary','normal','complex','symmetric', ...
    'parameter-dependent','eigenvalues','scalable'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite','<=',256});
if nargin<2||isempty(tau),tau=0;end
if nargin<3||isempty(gamma),gamma=1;end
if nargin<4||isempty(a),a=0;end
if nargin<5,b=[];end
if nargin<6||isempty(output),output='matrix';end
validateattributes(tau,{'numeric'},{'scalar','real','finite'});
validateattributes(gamma,{'numeric'},{'scalar','real','positive','finite'});
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('toda_cayley:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','operator','spectral','jacobi','parameters'}))
    error('toda_cayley:Mode','Unsupported output mode.');
end
if isstruct(a)
    if ~isempty(b)||~isscalar(a)||~isfield(a,'diagonal') ...
            ||~isfield(a,'subdiagonal')
        error('toda_cayley:Recurrence', ...
            'Supply a struct with diagonal/subdiagonal and empty B.');
    end
    d=a.diagonal;
    s=a.subdiagonal;
    if ~isnumeric(d)||~isreal(d)||numel(d)~=n ...
            ||any(~isfinite(d(:)))||~isvector(d) ...
            ||~isnumeric(s)||~isreal(s)||numel(s)~=n-1 ...
            ||any(~isfinite(s(:)))||any(s(:)<=0) ...
            ||(~isvector(s)&&~isempty(s))
        error('toda_cayley:Recurrence', ...
            'A Jacobi recurrence needs real diagonal and positive subdiagonal.');
    end
    d=d(:);s=s(:);
else
    if isempty(b),b=0;end
    [d,s]=js_jacobi_recurrence(a,b,n);
end
J0=diag(d);
if n>1
    J0=J0+diag(s,1)+diag(s,-1);
end
% eig(J0) gives independently fixed spectral data, even at large TAU.
[V,X]=eig(J0);
x=diag(X);
w0=abs(V(1,:).').^2;
w0=w0/sum(w0);
logw=log(w0)+2*tau*x;
logw=logw-max(logw);
w=exp(logw);
w=w/sum(w);
z=(x-1i*gamma)./(x+1i*gamma);
if strcmp(output,'spectral')
    U=unitary_spectral_reference('toda_cayley',n);
    U.finite.kind='exact';
    U.finite.eigenvalues=@()z;
    U.finite.moments=@(k)moments(z,ones(n,1)/n,k);
    U.finite.atoms=struct('angles',angle(z), ...
        'weights',ones(n,1)/n);
    U.cyclic=struct('kind','exact','angles',angle(z), ...
        'weights',w,'moments',@(k)moments(z,w,k));
    U.toda=struct('flow_time',tau,'cayley_scale',gamma, ...
        'jacobi_eigenvalues',x,'initial_weights',w0, ...
        'cyclic_weights',w,'weights_underflowed',any(w==0));
    U.reference='Symes, Physica D 4 (1982); Aurentz et al., ETNA 46 (2017)';
    U.notes=['Unitary symmetric Cayley image of finite Jacobi Toda QR flow. ' ...
        'Finite trace spectrum fixed; e1 cyclic weights deform. ' ...
        'At extreme times weights can underflow in floating point.'];
    return
end
% exp(tau*J0) and exp(tau*(J0-c*I)) have the SAME Q factor when
% c is real; using c=max(tau*x)/tau avoids exponential overflow.
if tau==0
    E=eye(n);
else
    shift=max(tau*x);
    E=expm(tau*J0-shift*eye(n));
end
[Q,R]=qr(E);
r=diag(R);
if any(abs(r)<=eps(max(1,norm(E,'fro')))*n)
    error('toda_cayley:QRBreakdown', ...
        'QR exponential is numerically rank deficient; reduce |TAU|.');
end
Q=Q*diag(r./abs(r));  % enforce positive diagonal R
M=Q'*J0*Q;
dnew=real(diag(M));
if n>1
    snew=real((diag(M,1)+diag(M,-1))/2);
    if any(snew<=0)
        error('toda_cayley:QRBreakdown', ...
            'Numerical Jacobi subdiagonal lost positivity.');
    end
    J=spdiags(dnew,0,n,n) ...
        +sparse(1:n-1,2:n,snew,n,n) ...
        +sparse(2:n,1:n-1,snew,n,n);
else
    J=sparse(dnew);
end
if strcmp(output,'jacobi'),U=J;return;end
if strcmp(output,'parameters')
    U=struct('initial_jacobi',sparse(J0),'jacobi',J, ...
        'flow_time',tau,'cayley_scale',gamma,'nodes',z, ...
        'jacobi_eigenvalues',x,'initial_weights',w0, ...
        'weights',w,'orthogonal_factor',Q, ...
        'tridiagonal_drift',norm(M-full(J),'fro'));
elseif strcmp(output,'operator')
    plusFactor=decomposition(J+1i*gamma*speye(n),'lu');
    minusFactor=decomposition(J-1i*gamma*speye(n),'lu');
    U=@(X,flag)um_cayley_apply(X,flag,plusFactor,minusFactor,gamma,n);
else
    K=full(J);I=eye(n);
    U=(K-1i*gamma*I)/(K+1i*gamma*I);
end
end

function m=moments(z,w,k)
if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
        ||any(k(:)~=fix(k(:)))
    error('toda_cayley:Moments','Moment powers must be finite integers.');
end
m=zeros(size(k));
for j=1:numel(k)
    m(j)=sum(w.*z.^k(j));
end
end
