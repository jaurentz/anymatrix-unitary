function [U,properties]=walsh_hadamard(n,ordering,output)
%WALSH_HADAMARD Normalized real Hadamard/Walsh transforms (order 2^k).
%   WALSH_HADAMARD(N,'natural') gives symmetric Sylvester order.
%   WALSH_HADAMARD(N,'sequency') gives increasing zero-crossing order
%   via bit-reversed Gray codes; default is 'sequency'.
%   WALSH_HADAMARD(N,ORDERING,'operator') uses the fast butterflies.
%   Applications: spread-spectrum communication and image coding.
%   Reference: K. G. Beauchamp, Applications of Walsh and Related
%   Functions, Academic Press, 1984.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [55] https://openlibrary.org/works/OL2903378W/Applications_of_Walsh_and_related_functions_with_an_introduction_to_sequency_theory
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('walsh_hadamard',n,ordering);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
m=round(log2(n));
if 2^m~=n,error('walsh_hadamard:Order','N must be a power of two.');end
if nargin<2 || isempty(ordering),ordering='sequency';end
if nargin<3,output='matrix';end
if nargin==2 && any(strcmpi(char(ordering),{'operator','implicit'}))
    output=ordering;ordering='sequency';
end
if ~(ischar(ordering)||isstring(ordering)) || ...
        ~any(strcmpi(char(ordering),{'natural','hadamard','sequency'}))
    error('walsh_hadamard:Ordering','Use natural or sequency ordering.');
end
if any(strcmpi(char(ordering),{'natural','hadamard'}))
    perm=(1:n)';
else
    k=uint64((0:n-1)');
    gray=bitxor(k,bitshift(k,-1));
    p=zeros(n,1);
    for b=1:m
        p=p+double(bitget(gray,b))*2^(m-b);
    end
    perm=p+1;
end
if um_is_operator(output)
    U=@(X,flag)um_walsh_apply(X,flag,perm);
else
    U=um_walsh_apply(eye(n),'notransp',perm);
end
end
