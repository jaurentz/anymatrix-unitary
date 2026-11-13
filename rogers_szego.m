function [U,properties] = rogers_szego(n,q,representation,beta,twists,output)
%ROGERS_SZEGO Rogers-Szego wrapped-normal real OPUC CMV generator.
%   Returns a finite N-by-N unitary multiplication matrix using the
%   measure -> Verblunsky -> terminal phase -> core factorization pipeline.
%   Supported forms: cmv (default), hessenberg, inverse_hessenberg, twisted.
%   Interior coefficients must lie in the open disk; |BETA|=1 closes.
%   See documentation and the corresponding SZEGO_*_MEASURE function.
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(q),q=0.6;end
if nargin<3 || isempty(representation),representation='cmv';end
if nargin<4 || isempty(beta),beta=1;end
if nargin<5,twists=[];end
mu=szego_rogers_measure(q);
if nargin<6,output='matrix';end
U=szego_matrix(mu,n,representation,beta,twists,output);
end
