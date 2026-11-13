function [U,properties] = chalker_coddington(nx,ny,theta,phases)
%CHALKER_CODDINGTON Sparse finite Chalker-Coddington network on a torus.
%   U=CHALKER_CODDINGTON(NX,NY,THETA,PHASES) has order 4*NX*NY.
%   Directed link coordinates are (x,y) on a 2*NX by 2*NY periodic
%   grid, with MATLAB index 1+x+(2*NX)*y (0-based x,y).
%
%   Two sets of 2-in/2-out scatterers implement the standard oriented
%   checkerboard network topology (Joye, Rev. Math. Phys. 2026, eq. 3.3).
%   The 'even' and 'odd' scattering matrices in the stated port order are
%      S_even=[cos(theta), sin(theta); -sin(theta), cos(theta)]
%      S_odd =[sin(theta), cos(theta); -cos(theta), sin(theta)].
%   THETA may be scalar or [theta_even,theta_odd]; default pi/4.
%   PHASES defaults to zeros; it can be a real scalar or 4*NX*NY
%   real vector, assigning one phase to each input directed link.
%   For independent random link disorder, supply 2*pi*rand(4*NX*NY,1)
%   after setting the caller's desired RNG state.
%
%   U is sparse with at most two nonzeros per row and column. Each
%   unitary scattering map sends disjoint input channel pairs to disjoint
%   output channel pairs, followed by diagonal input-link phases.
%   At THETA=0 or pi/2 it is a phase-weighted permutation.
%   The directed bipartite lattice obeys C*U*C=-U, where
%   C(x,y)=(-1)^(x+y). Consequently the finite spectrum is pi-paired
%   (lambda and -lambda), and all odd trace moments vanish exactly.
%
%   References: Chalker and Coddington, J. Phys. C 21 (1988), 2665;
%   Ho and Chalker, Phys. Rev. B 54 (1996), 8708;
%   Joye, Rev. Math. Phys. (2026), 2650005, eq. (3.3).
properties={'square','unitary','normal','complex','sparse','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(nx,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(ny),ny=nx;end
if nargin<3 || isempty(theta),theta=pi/4;end
validateattributes(ny,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(theta,{'numeric'},{'real','finite','vector'});
if ~(numel(theta)==1 || numel(theta)==2)
    error('chalker_coddington:Angles','THETA must have one or two entries.');
end
if isscalar(theta),theta=[theta;theta];else,theta=theta(:);end
N=4*nx*ny;
if nargin<4 || isempty(phases),phases=zeros(N,1);end
validateattributes(phases,{'numeric'},{'vector','real','finite'});
if isscalar(phases),phases=repmat(phases,N,1);end
if numel(phases)~=N
    error('chalker_coddington:Phases', ...
          'PHASES must be scalar or have 4*NX*NY entries.');
end
phases=phases(:);
A=[cos(theta(1)),sin(theta(1));-sin(theta(1)),cos(theta(1))];
B=[sin(theta(2)),cos(theta(2));-cos(theta(2)),sin(theta(2))];
X=2*nx;Y=2*ny;
id=@(x,y) 1+mod(x,X)+X*mod(y,Y);
row=zeros(2*N,1);
col=zeros(2*N,1);
val=zeros(2*N,1);
pos=0;
for j=0:nx-1
    for k=0:ny-1
        % Even node S_(2j,2k): columns (2j,2k),(2j+1,2k-1);
        % rows (2j+1,2k),(2j,2k-1), in the port order of eq. (3.3).
        c1=id(2*j,2*k);
        c2=id(2*j+1,2*k-1);
        r1=id(2*j+1,2*k);
        r2=id(2*j,2*k-1);
        q=pos+(1:4);
        row(q)=[r1;r2;r1;r2];
        col(q)=[c1;c1;c2;c2];
        val(q)=A(:);
        pos=pos+4;
        % Odd node S_(2j+1,2k): columns (2j+2,2k+1),(2j+1,2k);
        % rows (2j+1,2k+1),(2j+2,2k).
        c1=id(2*j+2,2*k+1);
        c2=id(2*j+1,2*k);
        r1=id(2*j+1,2*k+1);
        r2=id(2*j+2,2*k);
        q=pos+(1:4);
        row(q)=[r1;r2;r1;r2];
        col(q)=[c1;c1;c2;c2];
        val(q)=B(:);
        pos=pos+4;
    end
end
U=sparse(row,col,val.*exp(1i*phases(col)),N,N);
end
