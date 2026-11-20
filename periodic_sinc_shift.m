function [U,properties]=periodic_sinc_shift(n,delta,output)
%PERIODIC_SINC_SHIFT Real orthogonal periodic fractional-grid shift.
%   U=PERIODIC_SINC_SHIFT(N,DELTA) is the N-point trigonometric
%   interpolation/translation matrix, N positive and ODD.
%   DELTA is the shift in grid samples (default 0.375).
%   With centered integer frequencies k=-(N-1)/2:((N-1)/2),
%       U = F'*diag(exp(-2*pi*1i*k*DELTA/N))*F,
%   using the negative-exponent unitary DFT F.
%   This is a real orthogonal circulant matrix for every real
%   DELTA, and its eigenvalues are exactly the listed phases.
%   U(i,j) is a translated periodic Dirichlet cardinal kernel,
%   a finite periodic analogue of the sinc interpolation transform.
%   For DELTA=1 it is the forward cyclic shift.
%   U=PERIODIC_SINC_SHIFT(N,DELTA,'operator') returns a function
%   handle supporting flags 'notransp','transp','ctransp','size'.
%   That implementation uses FFTs in O(N log N) per right-hand side.
%   Even N is intentionally excluded: real orthogonal fractional
%   shifts require an explicit Nyquist-frequency convention.
properties={'square','unitary','normal','orthogonal','real', ...
    'circulant','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if mod(n,2)==0
    error('periodic_sinc_shift:Order','N must be odd.');
end
if nargin<2||isempty(delta),delta=.375;end
if nargin<3||isempty(output),output='matrix';end
if ischar(delta)||isstring(delta)
    if nargin>2,error('periodic_sinc_shift:Mode','Ambiguous output mode.');end
    output=delta;delta=.375;
end
validateattributes(delta,{'numeric'},{'scalar','real','finite'});
k=(0:n-1)';
k(k>(n-1)/2)=k(k>(n-1)/2)-n;
phase=exp(-2*pi*1i*k*delta/n);
if ~(ischar(output)||isstring(output))
    error('periodic_sinc_shift:Mode','Output must be matrix or operator.');
end
switch lower(char(output))
    case 'operator'
        U=@(X,flag)um_periodic_sinc_apply(X,flag,phase,n);
    case 'matrix'
        U=real(ifft(bsxfun(@times,phase,fft(eye(n),[],1)),[],1));
    otherwise
        error('periodic_sinc_shift:Mode','Output must be matrix or operator.');
end
end
