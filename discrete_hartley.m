function [U,properties] = discrete_hartley(n,output)
%DISCRETE_HARTLEY Real orthonormal discrete Hartley transform.
%   U(j+1,k+1) = [cos(2*pi*j*k/n)+sin(2*pi*j*k/n)]/sqrt(n).
%   U'=U and U^2=I, hence the only eigenvalues are +1 and -1.
%   Reference: R. N. Bracewell, JOSA 73 (1983), 1832-1835.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [47] https://doi.org/10.1364/JOSA.73.001832
properties={'square','unitary','normal','orthogonal','real', ...
    'symmetric','involutory','scalable','eigenvalues'};
% Optional analytic spectral-reference mode.
if nargin>=2 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('discrete_hartley',n);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2,output='matrix';end
if um_is_operator(output)
    U=@(X,flag)um_hartley_apply(X,flag,n);
    return
end
[J,K]=ndgrid(0:n-1,0:n-1);
a=2*pi*(J.*K)/n;
U=(cos(a)+sin(a))/sqrt(n);
end
