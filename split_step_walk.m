function [U, properties] = split_step_walk(n,theta1,theta2,output)
%SPLIT_STEP_WALK Periodic split-step discrete-time quantum walk.
%   U = SPLIT_STEP_WALK(N,THETA1,THETA2) returns a sparse 2N-by-2N
%   matrix acting on a periodic chain with a two-state coin. The basis
%   is ordered (x=0,up), (x=0,down), (x=1,up), (x=1,down), ... .
%
%   U = S_minus*C(THETA2)*S_plus*C(THETA1), where
%   C(THETA) acts at each site as the real rotation
%                [cos(THETA) -sin(THETA)]
%                [sin(THETA)  cos(THETA)].
%   S_plus shifts up states one site right; S_minus shifts down states
%   one site left. Both shifts have periodic boundary conditions.
%   Default angles are THETA1 = 0.3 and THETA2 = 0.7 radians.
%
%   This convention produces an exactly real orthogonal matrix in
%   exact arithmetic; it is a special subclass of unitary walks.
%
%   Related reference (note that conventions differ in the literature):
%   T. Kitagawa, M. S. Rudner, E. Berg, and E. Demler, Exploring
%   topological phases with quantum walks, Phys. Rev. A 82 (2010),
%   033429. https://doi.org/10.1103/PhysRevA.82.033429
properties = {'unitary', 'orthogonal', 'real', 'square', ...
              'sparse', 'scalable', 'parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=4
    if (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
        U=unitary_spectral_reference('split_step_walk',n,theta1,theta2);
        return
    end
    error('split_step_walk:Mode','Only terminal spectral mode is supported for this generator.');
end
if nargin == 0
    U = [];
    return
end
validateattributes(n, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
if nargin < 2
    theta1 = 0.3;
end
if nargin < 3
    theta2 = 0.7;
end
validateattributes(theta1, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(theta2, {'numeric'}, {'scalar', 'real', 'finite'});
col = 1:2*n;
row_plus = col;
row_minus = col;
row_plus(1:2:end) = 2*mod((0:n-1)+1,n)+1;
row_minus(2:2:end) = 2*mod((0:n-1)-1,n)+2;
S_plus = sparse(row_plus,col,1,2*n,2*n);
S_minus = sparse(row_minus,col,1,2*n,2*n);
R1 = [cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
R2 = [cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
C1 = kron(speye(n),sparse(R1));
C2 = kron(speye(n),sparse(R2));
U = S_minus*C2*S_plus*C1;
end
