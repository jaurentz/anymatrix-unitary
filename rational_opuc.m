function [U,properties] = rational_opuc(n,poles,schur,beta,output)
%RATIONAL_OPUC Finite unitary rational OPUC multiplication model.
%   U=RATIONAL_OPUC(N,POLES,SCHUR,BETA) applies the normal-parameter
%   operator Mobius transform to a finite unitary Hessenberg core matrix.
%   POLES contains N-1 prescribed exterior poles (abs(p)>1 or Inf).
%   SCHUR contains N-1 strict-disk recurrence parameters; BETA is a
%   terminal unimodular closure.  For POLES=Inf(N-1,1), U reduces to
%   the polynomial finite unitary Hessenberg construction.
%
%   Let a=[0;1./conj(POLES)], with 1/Inf defined as zero, A=diag(a),
%   D=(I-A'*A)^(1/2), and H the unitary core matrix.  Then
%       U = A + D*H*((I+A'*H)\D)
%   is exactly unitary for all valid parameters (operator Mobius map).
%   This is a finite rational-OPUC recurrence model. The poles and
%   recurrence coefficients jointly specify the associated measure;
%   it is NOT an arbitrary finite truncation of a rational Krylov pencil.
%   Finite eigenphase counting laws are generally not known in closed form.
%
%   References:
%   B. Simon, Orthogonal Polynomials on the Unit Circle, Parts 1-2.
%     https://doi.org/10.1090/coll054.1
%     https://doi.org/10.1090/coll/054.2
%   L. Velazquez, J. Funct. Anal. 254 (2008), 954-986,
%     https://doi.org/10.1016/j.jfa.2007.11.004
%   Bultheel, Cruz-Barroso, Lasarow, SIGMA 13 (2017), 090,
%     https://doi.org/10.3842/SIGMA.2017.090
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(poles),poles=inf(n-1,1);end
if nargin<3||isempty(schur),schur=zeros(n-1,1);end
if nargin<4||isempty(beta),beta=1;end
if nargin<5||isempty(output),output='matrix';end
if ~isnumeric(poles)||~isvector(poles)||numel(poles)~=n-1 ...
        ||any(isnan(poles(:)))||any(~isinf(poles(:)) & abs(poles(:))<=1)
    error('rational_opuc:Poles','POLES must be N-1 exterior poles (abs(p)>1 or Inf).');
end
if ~isnumeric(schur)||~isvector(schur)||numel(schur)~=n-1 ...
        ||any(~isfinite(schur(:)))||any(abs(schur(:))>=1)
    error('rational_opuc:Schur','SCHUR must be N-1 strict-disk coefficients.');
end
if ~isnumeric(beta)||~isscalar(beta)||~isfinite(beta) ...
        ||abs(abs(beta)-1)>1e-12
    error('rational_opuc:Boundary','BETA must have unit modulus.');
end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('rational_opuc:Mode','Output mode must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('rational_opuc',n);
    U.reference='Velazquez, J. Funct. Anal. 254 (2008), 954-986';
    U.notes=['The finite rational-OPUC Mobius model is unitary; ' ...
        'its generic eigenphase counting law is unknown. ' ...
        'Poles and recurrence parameters must not be confused with a ' ...
        'fixed underlying polynomial OPUC measure.'];
    return
end
if ~any(strcmp(output,{'matrix','explicit'}))
    error('rational_opuc:Mode','Only matrix and spectral modes are supported.');
end
a=zeros(n,1);
isfinitepole=~isinf(poles(:));
tmp=poles(:);
a(1+find(isfinitepole))=1./conj(tmp(isfinitepole));
A=diag(a);
D=diag(sqrt(1-abs(a).^2));
H=full(opuc_matrix(n,schur(:),beta,'hessenberg'));
U=A+D*H*((eye(n)+A'*H)\D);
end
