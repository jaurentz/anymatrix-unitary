function [U,properties] = time_frequency_shift(n,p,q,output)
%TIME_FREQUENCY_SHIFT Finite Weyl-Heisenberg modulation and translation.
%   U=M^P*T^Q, with T*e_j=e_(j+1 mod N) and
%   M*e_j=exp(2*pi*1i*j/N)*e_j, indexed from j=0.
%   This sparse monomial unitary admits exact cycle-based eigenvalues.
%   Application: discrete time-frequency shifts and Gabor analysis.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [51] https://doi.org/10.1109/78.564174
properties={'square','unitary','normal','sparse','complex', ...
    'scalable','eigenvalues','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=4 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('time_frequency_shift',n,p,q);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(p),p=1;end
if nargin<3 || isempty(q),q=1;end
if nargin<4,output='matrix';end
validateattributes(p,{'numeric'},{'scalar','integer','finite'});
validateattributes(q,{'numeric'},{'scalar','integer','finite'});
j=0:n-1;
row0=mod(j+q,n);
phase=exp(2*pi*1i*mod(p,n)*row0/n);
if um_is_operator(output)
    U=@(X,flag)um_monomial_apply(X,flag,row0+1,phase);
else
    U=sparse(row0+1,j+1,phase,n,n);
end
end
