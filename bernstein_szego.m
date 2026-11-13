function [U,properties] = bernstein_szego(n,a,representation,beta,twists,output)
%BERNSTEIN_SZEGO Finite OPUC matrix for the Poisson/Szego circle weight.
%   U=BERNSTEIN_SZEGO(N,A,REP,BETA,TWISTS) forms a finite multiplication
%   matrix from the probability measure with density
%     (1-|A|^2)/(2*pi*|exp(i*theta)-conj(A)|^2),  |A|<1.
%   This has Verblunsky coefficients alpha_0=A, alpha_k=0 (k>=1).
%   Default A=0.4, REP='cmv', closure BETA=1.
%   All the historic/twisted matrix representations share the OPUC
%   coefficient machinery; see JACOBI_SZEGO and SZEGO_MATRIX.
%   References: B. Simon, OPUC, AMS, 2005; DLMF section 18.33.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [5] https://doi.org/10.1090/coll054.1
%   [5] https://doi.org/10.1090/coll/054.2
%   [12] https://dlmf.nist.gov/18.33
properties={'square','unitary','normal','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
if nargin<2 || isempty(a),a=0.4;end
if nargin<3 || isempty(representation),representation='cmv';end
if nargin<4 || isempty(beta),beta=1;end
if nargin<5,twists=[];end
mu=js_bernstein_measure(a);
if nargin<6,output='matrix';end
U=js_measure_to_matrix(mu,n,representation,beta,twists,output);
end
