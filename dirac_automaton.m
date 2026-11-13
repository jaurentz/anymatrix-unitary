function [U,properties] = dirac_automaton(n,m,output)
%DIRAC_AUTOMATON Periodic one-dimensional Dirac quantum cellular automaton.
%   U=DIRAC_AUTOMATON(N,M) has size 2N by 2N and block ordering
%   [up on N sites; down on N sites]. With P the cyclic shift,
%       U=[nu*P, -1i*M*I; -1i*M*I, nu*P'], nu=sqrt(1-M^2).
%   0 <= M <= 1. Exact eigenvalues exp(+-1i*omega_k), where
%   cos(omega_k)=nu*cos(2*pi*k/N). The limiting eigenphase measure
%   has a mass-dependent band gap and square-root band edges (0<M<1).
%   U=DIRAC_AUTOMATON(N,M,'operator') gives a matrix-free action:
%       U(X,'notransp'), U(X,'transp'), U([],'size').
%   Reference: Bisio, D'Ariano, Tosini, Ann. Phys. 354 (2015) 244-264.
%   https://doi.org/10.1016/j.aop.2014.12.016
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [18] https://doi.org/10.1103/PhysRevA.93.012328
properties={'square','unitary','normal','sparse', ...
    'complex','eigenvalues','scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('dirac_automaton',n,m);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(m),m=0.4;end
validateattributes(m,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
if nargin<3,output='matrix';end
nu=sqrt((1-m)*(1+m));
if um_is_operator(output)
    U=@(X,flag)um_dirac_apply(X,flag,n,nu,m);
    return
end
P=cycle_shift(n);
I=speye(n);
U=[nu*P,-1i*m*I;-1i*m*I,nu*P'];
end
