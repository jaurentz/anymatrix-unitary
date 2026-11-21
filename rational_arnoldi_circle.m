function [U,properties]=rational_arnoldi_circle(n,poles,nodes,weights,output)
%RATIONAL_ARNOLDI_CIRCLE Finite unit-circle rational Arnoldi unitary model.
%   U=RATIONAL_ARNOLDI_CIRCLE(N,POLES,NODES,WEIGHTS,OUTPUT) constructs
%   orthonormal rational Krylov functions for a FINITE strictly positive
%   discrete probability measure sum_j w_j delta_(z_j) on the unit circle.
%   Distinct unimodular NODES z_j and positive WEIGHTS w_j are accepted.
%   POLES is a length N-1 sequence in the extended plane, allowing
%   finite poles inside or outside the unit disk and Inf (polynomial
%   step). Poles on or numerically extremely close to the circle are
%   rejected. Default nodes are N distinct equally spaced circle points.
%
%   The starting sample vector is q_1=sqrt(w). For j=1,...,N-1:
%     pole Inf: v=Z*q_j;
%     |pole|>1: v=Z*(I-Z/pole)^(-1)*q_j;
%     |pole|<1: v=(Z-pole I)^(-1)*q_j;
%   where Z=diag(z_j). The exterior expression is scaled to recover
%   the polynomial limit continuously when pole->Inf. Two-pass modified
%   Gram--Schmidt yields Q(:,1:j+1), with positive subdiagonal h_(j+1).
%   The final Q is square and unitary in exact arithmetic.
%
%   The N-by-N unitary finite CLOSURE returned by default is
%       U=Q'*Z*Q.
%   Thus spec(U)=NODES and the spectral weights of the first canonical
%   basis vector are exactly WEIGHTS, while trace weights remain 1/N.
%   This is a genuine exact finite weighted multiplication model,
%   not a generally nonunitary compression of an infinite operator.
%
%   The 'pencil' interface provides the **RECTANGULAR** N-by-(N-1)
%   upper-Hessenberg rational Arnoldi arrays H,K with U*K=H:
%      Inf: H(:,j)=h,        K(:,j)=e_j;
%      exterior p: H(:,j)=h, K(:,j)=e_j+h/p;
%      interior p: H(:,j)=e_j+p*h, K(:,j)=h.
%   Here h is the coefficient vector in v=Q*h. H,K are NOT themselves
%   square or unitary. (They are the rational recurrence, not an
%   arbitrary square generalized pencil.) The separate 'basis' output
%   contains Q, the sampled rational functions phi_j(z_i)=Q(i,j)/sqrt(w_i),
%   H,K and exact spectral data. 'spectral' returns the exact finite
%   uniform trace and distinct e1-cyclic moments/atoms.
%
%   OUTPUT: 'matrix' (default), 'basis', 'pencil', 'parameters', 'spectral'.
%   No 'operator': dense Q has O(N^2) storage and a handle would provide
%   no subquadratic benefit. The algorithm requires full (nonbreakdown)
%   rational Arnoldi; mixed poles can cause early exact breakdown or
%   near-breakdown and then the routine throws a diagnostic error.
%   N<=128. In contrast, RATIONAL_OPUC uses an operator Mobius transform
%   of a Schur-parameter Hessenberg unitary; these are different
%   parameterizations and not asserted to encode the same measure.
%
%   [162] Van Buggenhout, Van Barel, Vandebril, Numer. Algorithms 89
%   (2022), rational Arnoldi/Hessenberg-pencil inverse problems:
%   https://doi.org/10.1007/s11075-021-01125-6
%   [9] Bultheel, Cantero, Cruz-Barroso, unit-circle quadrature survey:
%   https://doi.org/10.1016/j.cam.2014.11.002
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',128});
if nargin<2||isempty(poles),poles=inf(n-1,1);end
if nargin<3||isempty(nodes)
    nodes=exp(2i*pi*((0:n-1)'+.17)/n);
end
if nargin<4||isempty(weights),weights=(1:n)';end
if nargin<5||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('rational_arnoldi_circle:Mode','OUTPUT must be scalar text.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','basis','pencil','parameters','spectral'}))
    error('rational_arnoldi_circle:Mode', ...
        'Use matrix, basis, pencil, parameters or spectral.');
end
if ~isnumeric(poles)||numel(poles)~=n-1|| ...
        (n>1 && ~isvector(poles)) || any(isnan(poles(:))) || ...
        any(~(isfinite(poles(:))|isinf(poles(:)))) || ...
        any(abs(abs(poles(isfinite(poles)))-1)<=1e-11)
    error('rational_arnoldi_circle:Poles', ...
        'POLES must be N-1 finite off-circle poles or Inf.');
end
poles=poles(:);
if ~isnumeric(nodes)||numel(nodes)~=n ||~isvector(nodes)|| ...
        any(~isfinite(nodes(:)))|| ...
        any(abs(abs(nodes(:))-1)>1e-11)
    error('rational_arnoldi_circle:Nodes', ...
        'NODES must contain N finite unit-circle nodes.');
end
nodes=nodes(:)./abs(nodes(:));
for j=2:n
    if min(abs(nodes(j)-nodes(1:j-1)))<=128*n*eps
        error('rational_arnoldi_circle:Distinct', ...
            'NODES must be distinct and numerically separated.');
    end
end
if ~isnumeric(weights)||~isreal(weights)||numel(weights)~=n|| ...
        ~isvector(weights)||any(~isfinite(weights(:)))|| ...
        any(weights(:)<=0)
    error('rational_arnoldi_circle:Weights', ...
        'WEIGHTS must contain N finite strictly positive reals.');
end
weights=weights(:);weights=weights/sum(weights);
if any(~isfinite(weights))||any(weights<=0)
    error('rational_arnoldi_circle:Weights', ...
        'Normalized positive weights must be representable.');
end
if strcmp(output,'spectral')
    U=unitary_spectral_reference('rational_arnoldi_circle',n);
    U.finite.kind='exact';
    U.finite.eigenvalues=@()nodes;
    U.finite.moments=@(k)weighted_moments(nodes,ones(n,1)/n,k);
    U.finite.atoms=struct('angles',angle(nodes), ...
        'weights',ones(n,1)/n);
    U.cyclic=struct('kind','exact','angles',angle(nodes), ...
        'weights',weights,'moments',@(k)weighted_moments(nodes,weights,k));
    U.poles=poles;
    U.reference='Van Buggenhout--Van Barel--Vandebril, Numerical Algorithms (2022)';
    U.notes=['A FINITE unitary closure Q''*diag(nodes)*Q from rational ', ...
        'Arnoldi weighted samples. Trace weights 1/N are different ', ...
        'from the exact e1-cyclic probability weights. The rational ', ...
        'recurrence H,K is a rectangular pencil, not a unitary matrix.'];
    return
end
Q=complex(zeros(n,n));
Q(:,1)=sqrt(weights);
H=complex(zeros(n,max(0,n-1)));
K=complex(zeros(n,max(0,n-1)));
stepnorm=zeros(max(0,n-1),1);
for j=1:n-1
    pole=poles(j);
    if isinf(pole)
        v=nodes.*Q(:,j);
        kind=0;
    elseif abs(pole)>1
        v=nodes./(1-nodes/pole).*Q(:,j);
        kind=1;
    else
        v=Q(:,j)./(nodes-pole);
        kind=2;
    end
    initial=norm(v);
    h=complex(zeros(n,1));
    for pass=1:2
        c=Q(:,1:j)'*v;
        h(1:j)=h(1:j)+c;
        v=v-Q(:,1:j)*c;
    end
    rho=norm(v);stepnorm(j)=rho;
    if ~isfinite(rho)||rho<=1e-12*initial
        error('rational_arnoldi_circle:Breakdown', ...
            'Rational Krylov basis broke down at step %d.',j);
    end
    Q(:,j+1)=v/rho;
    h(j+1)=rho;
    switch kind
        case 0
            K(j,j)=1;
            H(:,j)=h;
        case 1
            K(j,j)=1;
            K(:,j)=K(:,j)+h/pole;
            H(:,j)=h;
        case 2
            K(:,j)=h;
            H(j,j)=1;
            H(:,j)=H(:,j)+pole*h;
    end
end
T=Q'*(nodes.*Q);
if strcmp(output,'matrix'),U=T;return;end
recurrence=norm(T*K-H,'fro');
U=struct('dimension',n,'poles',poles,'nodes',nodes,'weights',weights, ...
    'Q',Q,'H',H,'K',K,'unitary_closure',T, ...
    'recurrence_residual',recurrence, ...
    'weighted_rational_basis',Q./sqrt(weights), ...
    'orthogonality_residual',norm(Q'*Q-eye(n),'fro'), ...
    'step_norms',stepnorm, ...
    'representation','finite weighted rational Arnoldi on unit-circle nodes', ...
    'pencil_size',[n,max(0,n-1)], ...
    'pencil_is_square',false);
if strcmp(output,'pencil')
    U=struct('H',H,'K',K,'dimension',n, ...
        'poles',poles,'recurrence_residual',recurrence, ...
        'pencil_size',[n,max(0,n-1)], ...
        'note','Rectangular rational Arnoldi recurrence T*K=H; H,K are not unitary.');
end
end

function m=weighted_moments(z,w,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
shape=size(k);p=reshape(k,1,[]);
m=reshape(sum(bsxfun(@times,w(:),bsxfun(@power,z(:),p)),1),shape);
end
