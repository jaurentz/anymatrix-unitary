function [U,properties] = circular_jacobi(n,a,eta,representation,beta,twists,output)
%CIRCULAR_JACOBI Generalized circular Jacobi complex-parameter OPUC CMV generator.
%   Returns a finite N-by-N unitary multiplication matrix using the
%   measure -> Verblunsky -> terminal phase -> core factorization pipeline.
%   Supported forms: cmv (default), hessenberg, inverse_hessenberg, twisted.
%   Interior coefficients must lie in the open disk; |BETA|=1 closes.
%   See documentation and the corresponding SZEGO_*_MEASURE function.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [5] https://doi.org/10.1090/coll054.1
%   [5] https://doi.org/10.1090/coll/054.2
%   [12] https://dlmf.nist.gov/18.33
%   [27] https://doi.org/10.1093/imrn/rnp092
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(a),a=0.3;end
if nargin<3 || isempty(eta),eta=0.5;end
if nargin<4 || isempty(representation),representation='cmv';end
if nargin<5 || isempty(beta),beta=1;end
if nargin<6,twists=[];end
mu=szego_circular_jacobi_measure(a,eta);
if nargin<7,output='matrix';end
U=szego_matrix(mu,n,representation,beta,twists,output);
end
