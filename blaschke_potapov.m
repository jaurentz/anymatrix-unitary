function [U,properties] = blaschke_potapov(n,poles,directions,omega,U0,output)
%BLASCHKE_POTAPOV Finite rational-inner (paraunitary) matrix response.
%   U=BLASCHKE_POTAPOV(N,POLES,DIRECTIONS,OMEGA,U0) evaluates
%      U0 * prod_j [I+(b_j(z)-1)*v_j*v_j'], z=exp(1i*OMEGA),
%      b_j(z)=(z-POLES(j))/(1-conj(POLES(j))*z), |POLES(j)|<1.
%   DIRECTIONS is N-by-M, M=numel(POLES); columns are normalized.
%   U0 is unitary (default I). Default factors are noncommuting for N>1.
%   Every result is N-by-N unitary at real OMEGA. A pointwise frequency
%   response is NOT automatically a finite sampled filter-bank operator.
%   OUTPUT='operator' provides the exact low-rank product action.
%
%   D. Alpay, P. Jorgensen, I. Lewkowicz, Opuscula Math. 36 (2016).
%   https://doi.org/10.7494/OpMath.2016.36.6.695
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(poles),poles=[.15+.2i;-.3+.1i];end
if ~isnumeric(poles)||~isvector(poles)||any(~isfinite(poles(:))) ...
        ||any(abs(poles(:))>=1)
    error('blaschke_potapov:Poles','POLES must lie strictly inside the unit disk.');
end
poles=poles(:);m=numel(poles);
if nargin<3||isempty(directions)
    directions=zeros(n,m);
    for j=1:m
        directions(mod(j-1,n)+1,j)=1;
        directions(mod(j,n)+1,j)=directions(mod(j,n)+1,j)+.5i;
    end
end
if ~isnumeric(directions)||~isequal(size(directions),[n m]) ...
        ||any(~isfinite(directions(:)))
    error('blaschke_potapov:Directions','DIRECTIONS must be finite N-by-M.');
end
for j=1:m
    nv=norm(directions(:,j));
    if nv==0,error('blaschke_potapov:Directions','No direction may vanish.');end
    directions(:,j)=directions(:,j)/nv;
end
if nargin<4||isempty(omega),omega=0;end
validateattributes(omega,{'numeric'},{'scalar','real','finite'});
if nargin<5||isempty(U0),U0=eye(n);end
if ~isnumeric(U0)||~isequal(size(U0),[n n]) ...
        ||any(~isfinite(U0(:)))||norm(U0'*U0-eye(n),'fro')>5e-11*max(1,n)
    error('blaschke_potapov:U0','U0 must be an N-by-N unitary matrix.');
end
if nargin<6||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('blaschke_potapov:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('blaschke_potapov',n);
    U.reference='Alpay, Jorgensen and Lewkowicz, Opuscula Math. 36 (2016)';
    U.notes=['Exactly unitary pointwise on the unit circle; no general ' ...
        'finite eigenphase law is known. OMEGA is a frequency parameter.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('blaschke_potapov:Mode','Mode must be matrix, operator or spectral.');
end
z=exp(1i*omega);
b=(z-poles)./(1-conj(poles)*z);
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_blaschke_apply(X,flag,n,b,directions,U0);
    return
end
U=U0;
for j=1:m
    v=directions(:,j);
    U=U+(b(j)-1)*(U*v)*v';
end
end
