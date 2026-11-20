function [U,properties]=complex_hadamard(n,m,twists,output)
%COMPLEX_HADAMARD Normalized Dita-deformed Fourier tensor construction.
%   U=COMPLEX_HADAMARD(N,M,TWISTS) has order N=M*R, with M dividing N.
%   Put F_M and F_R for the normalized DFTs, and
%      T = blkdiag(D_1*F_R,...,D_M*F_R), D_j=diag(exp(1i*TWISTS(:,j))).
%   The unitary complex Hadamard matrix is U=kron(F_M,I_R)*T.
%   TWISTS is R-by-M real; every entry of U has modulus 1/sqrt(N).
%   Empty TWISTS use deterministic nontrivial phase deformation when
%   M,R>1. M defaults to the smallest proper divisor (or 1 if prime).
%   M=1 with all-zero twists recovers the normalized Fourier matrix.
%   MODE='operator' uses two FFT axes; forward and adjoint cost
%   O(N*(log M+log R)) without assembling the dense matrix.
%   Generic deformed eigenphases are not available analytically.
%
%   P. Dita, Some results on the parametrization of complex Hadamard
%   matrices, J. Phys. A 37 (2004), 5355-5374.
%   https://doi.org/10.1088/0305-4470/37/20/008
%   W. Tadej and K. Zyczkowski, Open Syst. Inf. Dyn. 13 (2006).
%   https://doi.org/10.1007/s11080-006-8220-2
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(m)
    m=1;
    for k=2:floor(sqrt(double(n)))
        if mod(n,k)==0,m=k;break;end
    end
end
validateattributes(m,{'numeric'},{'scalar','integer','positive','finite'});
if mod(n,m)~=0
    error('complex_hadamard:Factor','M must be a positive divisor of N.');
end
r=n/m;
if nargin<3||isempty(twists)
    twists=(2*pi/(n+1))*((0:r-1)'*(0:m-1));
end
if ~isnumeric(twists)||~isreal(twists)||~isequal(size(twists),[r,m]) ...
        ||any(~isfinite(twists(:)))
    error('complex_hadamard:Twists','TWISTS must be a finite real (N/M)-by-M matrix.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('complex_hadamard:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('complex_hadamard',n);
    U.reference='Dita, J. Phys. A 37 (2004), doi:10.1088/0305-4470/37/20/008';
    if m==1 && ~any(twists(:))
        U=discrete_fourier(n,'spectral');
        U.family='complex_hadamard';
        U.reference='Normalized Fourier special case of the Dita construction.';
    else
        U.notes=['All entries have modulus 1/sqrt(N), but the generic ' ...
            'finite eigenphase counting measure is not known analytically.'];
    end
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('complex_hadamard:Mode','Supported modes: matrix, operator, spectral.');
end
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_hadamard_apply(X,flag,n,m,twists);
    return
end
F_m=exp(-2i*pi*((0:m-1)'*(0:m-1))/m)/sqrt(m);
F_r=exp(-2i*pi*((0:r-1)'*(0:r-1))/r)/sqrt(r);
U=zeros(n,n);
for j=1:m
    T=bsxfun(@times,exp(1i*twists(:,j)),F_r);
    cols=(j-1)*r+(1:r);
    for i=1:m
        rows=(i-1)*r+(1:r);
        U(rows,cols)=F_m(i,j)*T;
    end
end
end
