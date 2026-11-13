function [U,properties] = jacobi_szego(n,a,b,representation,beta,twists,output)
%JACOBI_SZEGO Jacobi-Szego multiplication matrices via Verblunsky cores.
%   U=JACOBI_SZEGO(N,A,B,REPRESENTATION,BETA,TWISTS) yields an N-by-N
%   unitary matrix from the Szego lift of the Jacobi probability measure
%   (1-t)^A(1+t)^B dt on [-1,1], with A,B>-1.
%   The circle density (with respect to dtheta) is proportional to
%   (1-cos(theta))^(A+1/2) * (1+cos(theta))^(B+1/2).
%
%   Default A=B=0 (Legendre), REPRESENTATION='cmv', BETA=1.
%   Representations: 'cmv', 'hessenberg', 'inverse_hessenberg', 'twisted'.
%   For 'twisted', TWISTS is a binary length N-2 snake pattern.
%   Any BETA on the unit circle is permitted; BETA=+1 or -1 gives
%   real orthogonal matrices, for all A,B.
%
%   Examples:
%     U=jacobi_szego(20,0,0,'cmv',1);
%     H=jacobi_szego(20,0,0,'hessenberg',1);
%     S=jacobi_szego(20,0,0,'twisted',exp(1i*.3),[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]);
%
%   References: Simon (2005), OPUC on One Foot, Section 11;
%   Bultheel-Cantero-Cruz-Barroso (2015), JCAM 284, 78-100.
%   https://dlmf.nist.gov/18.33
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [5] https://doi.org/10.1090/coll054.1
%   [5] https://doi.org/10.1090/coll/054.2
%   [7] https://doi.org/10.1016/j.cam.2006.10.033
%   [9] https://doi.org/10.1016/j.cam.2014.11.002
%   [11] https://doi.org/10.1090/S0273-0979-05-01075-X
properties = {'square','scalable','unitary','normal','parameter-dependent'};
if nargin==0
    U=[];return;
end
if nargin<2 || isempty(a), a=0; end
if nargin<3 || isempty(b), b=0; end
if nargin<4 || isempty(representation), representation='cmv'; end
if nargin<5 || isempty(beta), beta=1; end
if nargin<6, twists=[]; end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
mu=js_jacobi_measure(a,b); % A measure-first canonical entry point.
if nargin<7,output='matrix';end
U=js_measure_to_matrix(mu,n,representation,beta,twists,output);
end
