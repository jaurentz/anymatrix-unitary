function [U,properties] = chalker_coddington(nx,ny,theta,phases,output)
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
%   OUTPUT='operator' returns an O(N)-storage function handle using
%   a shared disjoint two-channel scattering action, including its
%   conjugate transpose. Default OUTPUT='matrix' remains sparse.
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
if nargin<5,output='matrix';end
implicit=um_is_operator(output);
[inputs,outputs,gates]=um_cc_geometry(nx,ny,theta);
phase=exp(1i*phases);
if implicit
    U=@(X,flag)um_scattering_apply(X,flag,inputs,outputs,gates,phase);
    return
end
% Assemble exactly the same two-channel network as the implicit action.
rows=[outputs(1,:);outputs(2,:);outputs(1,:);outputs(2,:)];
cols=[inputs(1,:);inputs(1,:);inputs(2,:);inputs(2,:)];
values=gates.*phase(cols);
U=sparse(rows(:),cols(:),values(:),N,N);
end
