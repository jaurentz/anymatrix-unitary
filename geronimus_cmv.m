function [U,properties] = geronimus_cmv(n,a,representation,beta,twists,output)
%GERONIMUS_CMV Geronimus constant-coefficient real or complex OPUC CMV generator.
%   Returns a finite N-by-N unitary multiplication matrix using the
%   measure -> Verblunsky -> terminal phase -> core factorization pipeline.
%   Supported forms: cmv (default), hessenberg, inverse_hessenberg, twisted.
%   Interior coefficients must lie in the open disk; |BETA|=1 closes.
%   See documentation and the corresponding SZEGO_*_MEASURE function.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [5] https://doi.org/10.1090/coll054.1
%   [5] https://doi.org/10.1090/coll/054.2
%   [7] https://doi.org/10.1016/j.cam.2006.10.033
properties={'square','unitary','normal','scalable','parameter-dependent'};
% Optional spectral-reference mode; preserve legacy output properties.
if nargin>=6 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('geronimus_cmv',n,a,representation,beta,twists);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(a),a=0.3+0.2i;end
if nargin<3 || isempty(representation),representation='cmv';end
if nargin<4 || isempty(beta),beta=1;end
if nargin<5,twists=[];end
mu=szego_geronimus_measure(a);
if nargin<6,output='matrix';end
U=szego_matrix(mu,n,representation,beta,twists,output);
end
