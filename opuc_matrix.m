function [U,properties] = opuc_matrix(n,alpha,beta,representation,twists,output)
%OPUC_MATRIX General finite unitary matrix from Verblunsky coefficients.
%   U=OPUC_MATRIX(N,ALPHA,BETA,REPRESENTATION,TWISTS) creates an N-by-N
%   unitary matrix using N-1 interior Verblunsky coefficients, a terminal
%   unit-circle beta, and a core factorization scheme.
%   N=1 is supported. Default representation is 'cmv'.
%   Options: 'hessenberg', 'inverse_hessenberg', 'cmv', 'twisted'.
%   See JS_CORE_PRODUCT for the TWISTS convention.
%   No-argument call returns the Anymatrix property list.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [8] https://etna.ricam.oeaw.ac.at/volumes/2011-2020/vol44/abstract.php?pages=327-341
%   [9] https://doi.org/10.1016/j.cam.2014.11.002
%   [10] https://doi.org/10.1016/0377-0427(86)90169-X
properties = {'square','scalable','unitary','normal','parameter-dependent'};
% Optional spectral-reference mode; preserve legacy output properties.
if nargin>=6 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('opuc_matrix',n,alpha,beta,representation,twists);
    return
end
if nargin==0
    U=[];return;
end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(alpha), alpha=zeros(n-1,1); end
if nargin<3 || isempty(beta), beta=1; end
if nargin<4 || isempty(representation), representation='cmv'; end
if nargin<5, twists=[]; end
if numel(alpha)~=n-1
    error('opuc_matrix:AlphaSize','Need N-1 interior coefficients.');
end
if nargin<6,output='matrix';end
U=js_core_product(alpha,beta,representation,twists,output);
end
