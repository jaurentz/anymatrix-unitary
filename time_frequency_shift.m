function [U,properties] = time_frequency_shift(n,p,q)
%TIME_FREQUENCY_SHIFT Finite Weyl-Heisenberg modulation and translation.
%   U=M^P*T^Q, with T*e_j=e_(j+1 mod N) and
%   M*e_j=exp(2*pi*1i*j/N)*e_j, indexed from j=0.
%   This sparse monomial unitary admits exact cycle-based eigenvalues.
%   Application: discrete time-frequency shifts and Gabor analysis.
properties={'square','unitary','normal','sparse','complex', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(p),p=1;end
if nargin<3 || isempty(q),q=1;end
validateattributes(p,{'numeric'},{'scalar','integer','finite'});
validateattributes(q,{'numeric'},{'scalar','integer','finite'});
j=0:n-1;
row0=mod(j+q,n);
phase=exp(2*pi*1i*mod(p,n)*row0/n);
U=sparse(row0+1,j+1,phase,n,n);
end
