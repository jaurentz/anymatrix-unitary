function [U,properties,A,B] = mobius_unitary_pencil(n,a,mode)
%MOBIUS_UNITARY_PENCIL Strict unitary pencil from a disk automorphism.
%   U=MOBIUS_UNITARY_PENCIL(N,A) returns B\A with
%     A=P-A*I, B=I-conj(A)*P, P=CYCLE_SHIFT(N).
%   [U,PROPERTIES,ACOEF,BCOEF] exposes ACOEF-z*BCOEF;
%   ACOEF*ACOEF'=BCOEF*BCOEF', so inv(BCOEF)*ACOEF is unitary.
%   0<=abs(A)<1. Default A=.4. The exact eigenvalues are
%   (omega_j-A)/(1-conj(A)*omega_j), omega_j=exp(2i*pi*j/N).
%   The limiting equal-weight density is a Poisson kernel.
%   This pencil specialization is derived from the disk automorphism.
%   Background: Simon, OPUC, AMS 2005.
%   [5] https://doi.org/10.1090/coll054.1
%   [5] https://doi.org/10.1090/coll/054.2
properties={'square','unitary','normal','complex','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,U=[];A=[];B=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(a),a=.4;end
validateattributes(a,{'numeric'},{'scalar','finite'});
if abs(a)>=1,error('mobius_unitary_pencil:Parameter','Require abs(a)<1.');end
if nargin<3 || isempty(mode),mode='matrix';end
if strcmpi(char(mode),'spectral')
    U=ep_spectral_reference('mobius_unitary_pencil',n,a);A=[];B=[];return
end
if ~strcmpi(char(mode),'matrix')
    error('mobius_unitary_pencil:Mode','Supported modes: matrix, spectral.');
end
P=cycle_shift(n);I=speye(n);
A=P-a*I; B=I-conj(a)*P;
U=B\A;
end
