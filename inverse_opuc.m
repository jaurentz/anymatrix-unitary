function [U,properties]=inverse_opuc(n,nodes,weights,representation,output)
%INVERSE_OPUC Finite CMV/Hessenberg inverse spectral construction.
%   U=INVERSE_OPUC(N,NODES,WEIGHTS,REP,OUTPUT) returns a genuinely
%   unitary N-by-N CMV or unitary Hessenberg matrix with prescribed
%   DISTINCT eigenvalues NODES on the unit circle and strictly positive
%   e1 spectral weights WEIGHTS, internally normalized to sum one.
%   WEIGHTS are NOT the equal-weight trace eigenvalue distribution.
%   REP='cmv' (default), 'hessenberg' or 'inverse_hessenberg'.
%
%   Interior Verblunsky coefficients are recovered via finite weighted
%   Szego orthogonalization, normalizing the recursively evaluated
%   polynomial and reversed polynomial at each step. The terminal
%   phase is beta=(-1)^(N+1)*conj(prod(NODES)).
%   The output reuses JS_CORE_PRODUCT, never a nonunitary companion.
%   'parameters' returns alpha (N-1), beta and input nodes/weights.
%   'operator' provides exact core-product actions with adjoint.
%   'spectral' returns the EXACT finite trace moments of the given
%   nodes and a separate e1 weighted atomic measure in S.cyclic.
%   A repeated node, zero weight or numerically singular recurrence
%   is rejected; no claim of stability for arbitrarily clustered nodes.
%   [5] https://doi.org/10.1090/coll054.1
%   [5] https://doi.org/10.1090/coll/054.2
%   [9] https://doi.org/10.1016/j.cam.2014.11.002
properties={'square','unitary','normal','scalable','parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(nodes)
    nodes=exp(2i*pi*((0:n-1)'+.17)/n);
end
if ~isnumeric(nodes)||~isvector(nodes)||numel(nodes)~=n ...
        ||any(~isfinite(nodes(:))) ...
        ||any(abs(abs(nodes(:))-1)>1e-11)
    error('inverse_opuc:Nodes','NODES must be N unimodular numbers.');
end
nodes=nodes(:)./abs(nodes(:));
for j=2:n
    if min(abs(nodes(j)-nodes(1:j-1))) <= 128*eps*max(1,n)
        error('inverse_opuc:Distinct','Eigenphases must be distinct.');
    end
end
if nargin<3||isempty(weights),weights=ones(n,1)/n;end
if ~isnumeric(weights)||~isreal(weights)||~isvector(weights) ...
        ||numel(weights)~=n||any(~isfinite(weights(:))) ...
        ||any(weights(:)<=0)
    error('inverse_opuc:Weights','WEIGHTS must be N strictly positive reals.');
end
weights=weights(:);
weights=weights/sum(weights);
if any(weights==0)||any(~isfinite(weights))
    error('inverse_opuc:Weights','Positive weights must sum finitely.');
end
if nargin<4||isempty(representation),representation='cmv';end
if ~(ischar(representation)||(isstring(representation)&&isscalar(representation)))
    error('inverse_opuc:Representation','REP must be a scalar string.');
end
representation=lower(char(representation));
if ~any(strcmp(representation,{'cmv','hessenberg','inverse_hessenberg'}))
    error('inverse_opuc:Representation','Unsupported OPUC core representation.');
end
if nargin<5||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('inverse_opuc:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('inverse_opuc',n);
    U.finite.kind='exact';
    U.finite.eigenvalues=@()nodes;
    U.finite.moments=@(k)circle_moments(nodes,ones(n,1)/n,k);
    U.finite.atoms=struct('angles',angle(nodes),'weights',ones(n,1)/n);
    U.cyclic=struct('kind','exact','angles',angle(nodes), ...
        'weights',weights,'moments',@(k)circle_moments(nodes,weights,k));
    U.reference='Simon, OPUC (2005), and unit-circle matrix quadrature';
    U.notes=['Finite trace measure uses weights 1/N; supplied WEIGHTS ' ...
        'are the distinct e1 spectral measure, NOT the trace measure.'];
    return
end
alpha=inverse_verblunsky(nodes,weights);
beta=alpha(end);
switch output
    case 'parameters'
        U=struct('alpha',alpha(1:end-1),'beta',beta, ...
            'nodes',nodes,'weights',weights,'representation',representation, ...
            'dimension',n);
    case {'matrix','explicit'}
        U=js_core_product(alpha(1:end-1),beta,representation);
    case {'operator','implicit'}
        U=js_core_product(alpha(1:end-1),beta,representation,[],'operator');
    otherwise
        error('inverse_opuc:Mode','Use matrix, operator, spectral or parameters.');
end
end

function alpha=inverse_verblunsky(z,w)
% In exact arithmetic r_k=phi_k(z), s_k=phi_k^*(z), where phi_k
% is monic Szego polynomial scaled to unit norm; their common
% normalization does not change the reflection coefficient.
n=numel(z);r=ones(n,1);s=ones(n,1);
alpha=complex(zeros(n,1));
for k=1:n
    den=sum(w.*abs(s).^2);
    if den<=realmin||~isfinite(den)
        error('inverse_opuc:NearBreakdown','Ill-conditioned OPUC recurrence.');
    end
    c=sum(w.*conj(s).*z.*r)/den;   % c=conj(alpha_(k-1))
    if k==n
        beta=(-1)^(n+1)*conj(prod(z));
        if abs(abs(c)-1)>3e-6||abs(conj(c)-beta)>3e-6
            error('inverse_opuc:NearBreakdown', ...
                'Terminal phase inconsistent with finite support.');
        end
        alpha(k)=beta;
        return
    end
    if ~isfinite(c)||abs(c)>=1-1e-12
        error('inverse_opuc:NearBreakdown', ...
            'Interior Verblunsky coefficient too close to unit circle.');
    end
    alpha(k)=conj(c);
    nextR=z.*r-c.*s;
    nextS=s-conj(c).*z.*r;
    scale=sqrt((1-abs(c))*(1+abs(c)));
    r=nextR/scale;
    s=nextS/scale;
end
end

function m=circle_moments(z,w,k)
if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
        ||any(k(:)~=round(k(:)))
    error('inverse_opuc:Moments','Moment orders must be finite integers.');
end
m=zeros(size(k));
for j=1:numel(k)
    m(j)=sum(w.*z.^k(j));
end
end
