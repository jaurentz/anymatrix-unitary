function [U,properties]=periodic_sinc_shift(n,delta,output,variant)
%PERIODIC_SINC_SHIFT Periodic fractional-grid translation matrices.
%   Q=PERIODIC_SINC_SHIFT(N,DELTA) gives the REAL ORTHOGONAL
%   Dirichlet-cardinal translation for positive ODD N (unchanged API).
%   DELTA is measured in grid samples (default 0.375).
%
%   For EVEN N two explicit, mathematically different variants are
%   available with a fourth argument:
%     PERIODIC_SINC_SHIFT(N,DELTA,'matrix','complex')
%       Continuous complex UNITARY fractional shift. The unpaired
%       Nyquist Fourier frequency is +N/2, so Q(1) is the complete
%       forward cyclic permutation. Q(t) is generally nonreal.
%     PERIODIC_SINC_SHIFT(N,DELTA,'matrix','nyquist-frozen')
%       Continuous REAL ORTHOGONAL fractional shift with Nyquist
%       multiplier fixed at +1. Q(1) differs from the one-step cyclic
%       permutation; it is NOT its real orthogonal interpolation.
%
%   Each variant has Q(t)*Q(s)=Q(t+s), Q(0)=I, and Fourier eigenvalues
%   exp(-2*pi*1i*k*t/N), with the indicated frequency convention.
%   For even N, a continuous real orthogonal interpolation from I to
%   the one-step shift is impossible because det(P)=-1.
%
%   OUTPUT='operator' gives a function handle with 'notransp',
%   'transp'/'ctransp' (both conjugate adjoints), and 'size' flags.
%   OUTPUT='spectral' gives exact finite eigenvalues and moments.
%   Both matrix-free implementations use O(N log N) FFT products.
%   [1] https://openlibrary.org/books/OL4407364M/Circulant_matrices
%   [52] https://doi.org/10.1109/TSP.2008.925261
properties={'square','unitary','normal','orthogonal','real', ...
    'circulant','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(delta),delta=.375;end
if nargin<3||isempty(output),output='matrix';end
if ischar(delta)||isstring(delta)
    if nargin>2,error('periodic_sinc_shift:Mode','Ambiguous output mode.');end
    output=delta;delta=.375;
end
if nargin<4||isempty(variant),variant='real';end
validateattributes(delta,{'numeric'},{'scalar','real','finite'});
if ~(ischar(variant)||isstring(variant))
    error('periodic_sinc_shift:Variant','VARIANT must be a string.');
end
variant=lower(char(variant));
if ~any(strcmp(variant,{'real','complex','nyquist-frozen'}))
    error('periodic_sinc_shift:Variant','Unknown fractional shift variant.');
end
if mod(n,2)==0
    if strcmp(variant,'real')
        error('periodic_sinc_shift:Order', ...
            'Even N requires explicit complex or nyquist-frozen variant.');
    end
else
    if ~strcmp(variant,'real')
        error('periodic_sinc_shift:Variant', ...
            'The named complex and nyquist-frozen variants require even N.');
    end
end
if strcmp(variant,'complex')
    properties={'square','unitary','normal','complex', ...
        'circulant','scalable','parameter-dependent'};
end
if ~(ischar(output)||isstring(output))
    error('periodic_sinc_shift:Mode','Output must be a string.');
end
k=(0:n-1)';
if mod(n,2)==1
    k(k>(n-1)/2)=k(k>(n-1)/2)-n;
else
    k(k>n/2)=k(k>n/2)-n;
    if strcmp(variant,'nyquist-frozen')
        k(n/2+1)=0;
    end
end
phase=exp(-2*pi*1i*k*delta/n);
switch lower(char(output))
    case 'operator'
        if strcmp(variant,'complex')
            U=@(X,flag)um_fft_apply(X,flag,'circulant',phase,[]);
        else
            U=@(X,flag)um_periodic_sinc_apply(X,flag,phase,n);
        end
    case 'spectral'
        U=unitary_spectral_reference('periodic_sinc_shift',n,delta,variant);
    case 'matrix'
        U=ifft(bsxfun(@times,phase,fft(eye(n),[],1)),[],1);
        if ~strcmp(variant,'complex'),U=real(U);end
    otherwise
        error('periodic_sinc_shift:Mode', ...
            'Output must be matrix, operator or spectral.');
end
end
