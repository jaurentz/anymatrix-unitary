function [U,properties] = discrete_hartley(n)
%DISCRETE_HARTLEY Real orthonormal discrete Hartley transform.
%   U(j+1,k+1) = [cos(2*pi*j*k/n)+sin(2*pi*j*k/n)]/sqrt(n).
%   U'=U and U^2=I, hence the only eigenvalues are +1 and -1.
%   Reference: R. N. Bracewell, JOSA 73 (1983), 1832-1835.
properties={'square','unitary','normal','orthogonal','real', ...
    'symmetric','involutory','scalable','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
[J,K]=ndgrid(0:n-1,0:n-1);
a=2*pi*(J.*K)/n;
U=(cos(a)+sin(a))/sqrt(n);
end
