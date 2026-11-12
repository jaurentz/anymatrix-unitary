function [U,alpha] = szego_matrix(mu,n,representation,beta,twists)
%SZEGO_MATRIX Finite unitary multiplication matrix from a Szego measure.
%   [U,ALPHA] = SZEGO_MATRIX(MU,N,REP,BETA,TWISTS) takes the first N-1
%   Verblunsky coefficients of MU and closes at unit phase BETA.
%   REP can be 'cmv' (default), 'hessenberg', 'inverse_hessenberg', or
%   'twisted'. Optional TWISTS chooses the left/right insertion pattern.
%   The original measure and the finite e1 quadrature measure differ.
if nargin<3 || isempty(representation), representation='cmv'; end
if nargin<4 || isempty(beta), beta=1; end
if nargin<5, twists=[]; end
[U,alpha]=js_measure_to_matrix(mu,n,representation,beta,twists);
end
