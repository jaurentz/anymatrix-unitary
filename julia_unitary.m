function [U,properties] = julia_unitary(n,T,output)
%JULIA_UNITARY Halmos/Julia 2N-by-2N unitary from a contraction T.
%   U=[T, (I-T*T')^(1/2); (I-T'*T)^(1/2), -T'].
%   T must be an N-by-N contraction. Default T=0.5*CYCLE_SHIFT(N).
%   Defect square roots are computed by an SVD and are positive
%   semidefinite. This is a single finite block-unitary extension; it
%   does NOT imply compression identities for all powers of T.
%   Reference: B. Sz.-Nagy, Unitary Dilations of Hilbert Space Operators,
%   CBMS 19, AMS, 1974. https://doi.org/10.1090/CBMS/019
properties={'square','unitary','normal','scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('julia_unitary',n,T);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(T),T=0.5*full(cycle_shift(n));end
validateattributes(T,{'numeric'},{'size',[n n],'finite'});
[L,S,R]=svd(full(T));
s=diag(S);
if any(s>1+100*n*eps)
    error('julia_unitary:Contraction','The spectral norm of T must be <= 1.');
end
d=sqrt(max(0,1-min(s,1).^2));
DTstar=L*diag(d)*L';
DT=R*diag(d)*R';
U=[T,DTstar;DT,-T'];
end
