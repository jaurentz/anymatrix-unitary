function [U,properties]=time_skew_magnus(n,T,steps,frequency,amplitude,method,output)
%TIME_SKEW_MAGNUS Noncommuting time-dependent SO(n) integration benchmark.
%   U=TIME_SKEW_MAGNUS(N,T,STEPS,FREQUENCY,AMPLITUDE,METHOD,OUTPUT)
%   approximates Y'=A(t)*Y, Y(0)=I, for 0<=t<=T and
%       A(t)=A0+AMPLITUDE*sin(FREQUENCY*t)*A1.
%   A0 has skew-symmetric first offdiagonal links and A1 has
%   skew-symmetric second offdiagonal links. For N>=3, [A0,A1]~=0.
%   Both are real and sparse. The default method is 'cayley_midpoint'.
%   Every stage uses Aj=A((j-1/2)*h), h=T/STEPS.
%
%   'exp_midpoint': Qj=expm(h*Aj), second-order exponential midpoint
%       (the second-order commutator-free Magnus method).
%   'cayley_midpoint': Qj=(I-h*Aj/2)\(I+h*Aj/2), implicit midpoint
%       and a second-order Cayley approximation, NOT a 4th-order
%       modified Cayley--Magnus formula.
%   Ordered composition is U=Q_STEPS*...*Q_1.
%   Each stage is in SO(N) in exact arithmetic; the product is
%   exactly orthogonal and has det +1 (up to floating-point error).
%   Generic SO(2m) factors are NOT automatically symplectic.
%
%   OUTPUT='operator' applies stages and Hermitian adjoints, with
%       cached sparse LU Cayley solvers. The exp method caches DENSE
%       stage exponentials and does not claim a sparse expmv routine.
%   OUTPUT='data' gives A0,A1 and A(t) for independent verification.
%   OUTPUT='spectral' reports UNKNOWN nonautonomous trace eigenphases;
%       if AMPLITUDE=0, exact algebraic phases are derived from A0:
%       exp(T*eig(A0)) for exponential midpoint, and
%       ((1+h*eig(A0)/2)/(1-h*eig(A0)/2))^STEPS for Cayley.
%   N>=3 enforces a genuinely noncommuting available pair, except
%   AMPLITUDE=0 where the time-dependent contribution is disabled.
%
%   [122] Celledoni--Marthinsen--Owren (2003):
%   https://doi.org/10.1016/S0167-739X(02)00161-9
%   [101] Hairer--Lubich--Wanner (2006):
%   https://doi.org/10.1007/3-540-30666-8
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','>=',3,'finite'});
if nargin<2||isempty(T),T=1.25;end
if nargin<3||isempty(steps),steps=8;end
if nargin<4||isempty(frequency),frequency=2.3;end
if nargin<5||isempty(amplitude),amplitude=.7;end
if nargin<6||isempty(method),method='cayley_midpoint';end
if nargin<7||isempty(output),output='matrix';end
validateattributes(T,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(steps,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(frequency,{'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(amplitude,{'numeric'},{'scalar','real','finite'});
if ~(ischar(method)||(isstring(method)&&isscalar(method))) ...
        ||~any(strcmpi(char(method),{'exp_midpoint','cayley_midpoint'}))
    error('time_skew_magnus:Method','Choose exp_midpoint or cayley_midpoint.');
end
method=lower(char(method));
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('time_skew_magnus:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
j=(1:n-1)';
c0=.25+.03*cos(j);
A0=sparse(j,j+1,c0,n,n);A0=A0-A0';
j=(1:n-2)';
c1=.45+.025*sin(j);
A1=sparse(j,j+2,c1,n,n);A1=A1-A1';
h=T/steps;
if strcmp(output,'data')
    U=struct('A0',A0,'A1',A1, ...
        'at',@(time)A0+amplitude*sin(frequency*time)*A1, ...
        'midpoints',((1:steps)'-.5)*h,'step_size',h, ...
        'dimension',n,'method',method);
    return
end
if strcmp(output,'spectral')
    U=unitary_spectral_reference('time_skew_magnus',n);
    U.reference='Celledoni, Marthinsen and Owren (2003)';
    if amplitude==0
        lambda=eig(full(A0));
        if strcmp(method,'exp_midpoint')
            z=exp(T*lambda);
        else
            z=((1+.5*h*lambda)./(1-.5*h*lambda)).^steps;
        end
        z=z./abs(z);
        U.finite.kind='exact';
        U.finite.eigenvalues=@()z;
        U.finite.moments=@(k)um_skew_moments(z,k);
        U.finite.atoms=struct('angles',angle(z), ...
            'weights',ones(n,1)/n);
        U.notes=['Constant-generator spectral reduction using eig(A0) ' ...
            '(small numerical eigensolve); exact midpoint/Cayley ' ...
            'spectral mapping in exact arithmetic.'];
    else
        U.notes=['Time-dependent noncommuting skew generator: no ' ...
            'generic closed-form finite trace eigenphase law. ' ...
            'The integrator is not asserted to be standard-J symplectic.'];
    end
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('time_skew_magnus:Mode', ...
        'Choose matrix, operator, spectral, or data.');
end
if strcmp(method,'exp_midpoint') && n>128
    error('time_skew_magnus:Order', ...
        'Dense exponential midpoint limited to n<=128.');
end
if any(strcmp(output,{'matrix','explicit'})) && n>384
    error('time_skew_magnus:Order', ...
        'Explicit dense output limited to n<=384; use operator.');
end
stages=cell(steps,1);
I=speye(n);
for s=1:steps
    Aj=A0+amplitude*sin(frequency*(s-.5)*h)*A1;
    if strcmp(method,'exp_midpoint')
        stages{s}=expm(h*full(Aj));
    else
        C=.5*h*Aj;
        stages{s}=struct('minus',decomposition(I-C,'lu'), ...
            'plus',decomposition(I+C,'lu'));
    end
end
F=@(X,flag)um_skew_midpoint_apply(X,flag,n,stages,method);
if any(strcmp(output,{'operator','implicit'}))
    U=F;
else
    U=F(eye(n),'notransp');
    U=real(U);
end
end
