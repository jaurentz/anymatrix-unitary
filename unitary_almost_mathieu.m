function [W,properties] = unitary_almost_mathieu(n,lambda1,lambda2,flux,phase,mode)
%UNITARY_ALMOST_MATHIEU Finite periodic-window magnetic quantum walk.
%   W=UNITARY_ALMOST_MATHIEU(N,LAMBDA1,LAMBDA2,FLUX,PHASE)
%   constructs the sparse 2N-by-2N unitary S_lambda1*Q_lambda2,flux,phase
%   of Cedzich, Fillman and Ong (2023), equations (2.3)-(2.6).
%   Defaults: LAMBDA1=LAMBDA2=0.7, FLUX=(sqrt(5)-1)/2, PHASE=.13;
%   FLUX and PHASE are measured in cycles, not radians.
%   Q_j=[lambda2*cos(t_j)+i*r2, -lambda2*sin(t_j);
%        lambda2*sin(t_j),     lambda2*cos(t_j)-i*r2],
%   t_j=2*pi*(j*FLUX+PHASE), r2=sqrt(1-lambda2^2).
%   S_lambda maps |j,+> to lambda|j+1,+>+r1|j,-> and
%   |j,-> to lambda|j-1,->-r1|j,+>, r1=sqrt(1-lambda1^2).
%
%   For rational FLUX=p/q and q dividing N this is a finite periodic
%   cell. For irrational FLUX, the periodic boundary closes a finite
%   window; it is not an infinite irrational-periodic walk. The
%   subcritical/critical/supercritical spectral type statements are
%   properties of the infinite, not finite, operator.
%   'pencil' returns P.A=Q, P.B=S', P.U=W.
%   'spectral' returns unknown finite-trace spectral metadata.
%   https://doi.org/10.1007/s00220-023-04808-4
properties={'square','unitary','normal','complex','sparse','scalable', ...
            'parameter-dependent'};
if nargin==0,W=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','finite','positive'});
if nargin<2 || isempty(lambda1),lambda1=.7;end
if nargin<3 || isempty(lambda2),lambda2=.7;end
if nargin<4 || isempty(flux),flux=(sqrt(5)-1)/2;end
if nargin<5 || isempty(phase),phase=.13;end
validateattributes(lambda1,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
validateattributes(lambda2,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
validateattributes(flux,{'numeric'},{'scalar','real','finite'});
validateattributes(phase,{'numeric'},{'scalar','real','finite'});
if nargin<6 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('unitary_almost_mathieu:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if strcmp(mode,'spectral')
    W=unitary_spectral_reference('split_step_walk',n);
    W.family='unitary_almost_mathieu';
    W.reference='Cedzich, Fillman and Ong, Commun. Math. Phys. 403 (2023)';
    W.notes=['2N-state finite periodic window; infinite-volume regime ' ...
        'statements do not apply literally to this finite matrix.'];
    return
end
if ~any(strcmp(mode,{'matrix','pencil'}))
    error('unitary_almost_mathieu:Mode', ...
          'Supported modes: matrix, pencil, spectral.');
end
N=2*n;
j=(0:n-1)';a=2*j+1;b=a+1;
r1=sqrt(max(0,1-lambda1^2));
r2=sqrt(max(0,1-lambda2^2));
rows=[2*mod(j+1,n)+1;b;2*mod(j-1,n)+2;a];
cols=[a;a;b;b];
vals=[lambda1*ones(n,1);r1*ones(n,1); ...
      lambda1*ones(n,1);-r1*ones(n,1)];
S=sparse(rows,cols,vals,N,N);
t=2*pi*(j*flux+phase);c=lambda2*cos(t);s=lambda2*sin(t);
rows=[a;b;a;b];cols=[a;a;b;b];
vals=[c+1i*r2;s;-s;c-1i*r2];
Q=sparse(rows,cols,vals,N,N);
W=S*Q;
if strcmp(mode,'pencil')
    W=struct('A',Q,'B',S','U',W);
end
end
