function [U,properties]=abl_ladik_cmv(n,t,nodes,weights,harmonic,output)
%ABL_LADIK_CMV Exact finite defocusing Ablowitz--Ladik (Schur) Lax flow.
%   U=ABL_LADIK_CMV(N,T,NODES,WEIGHTS,HARMONIC,OUTPUT) constructs the
%   finite N-by-N unitary CMV matrix on the invariant isospectral leaf
%   with boundary Verblunsky coefficient ALPHA_(N-1)=-1. N must be EVEN
%   and at least 4 (the finite boundary reduction in Nenciu Section 4).
%   Z_j are distinct unimodular support nodes whose product yields
%   terminal alpha=-1; W_j>0 are initial e1 spectral weights.
%
%   We choose the Hamiltonian H_k = 2 Im Tr(C^k)/k, k=HARMONIC>=1.
%   The positive e1 spectral weights are evolved EXACTLY by
%     w_j(t)=w_j(0)*exp(-2*t*Re(z_j^k))/Z(t), Z=sum of numerator,
%   and the nodes stay fixed. The inverse finite OPUC construction
%   converts the weights to Verblunsky parameters and a five-diagonal
%   CMV matrix WITHOUT integrating their nonlinear differential ODE.
%
%   The corresponding finite Lax equation, with conventions of
%   Nenciu (2005), Theorem 3, Eq. (4.7), is
%     dC/dt=[C,B_k], B_k=(C^k)_+ - ((C^k)_+)^*,
%   where X_+ is strictly upper triangular X plus HALF its diagonal.
%   B_k is skew-Hermitian. 'lax' output gives C,B_k,[C,B_k], and the
%   predicted cyclic weight derivative. The k=1 real subclass obeys
%     d(alpha_j)/dt = (1-|alpha_j|^2)*(alpha_(j-1)-alpha_(j+1)),
%   with alpha_(-1)=alpha_(N-1)=-1. This is Schur/dmKdV flow, a
%   distinguished REAL Hamiltonian within the defocusing AL hierarchy;
%   it is not the usual complex NLS-like flow for Re Tr(C).
%
%   OUTPUT 'matrix' (default), 'operator', 'parameters', 'spectral',
%   'lax'. The operator is an O(N)-storage CMV core product after
%   O(N^2) inverse spectral setup; the 'lax' diagnostic forms dense
%   C^k for N<=128 only. Distinct nodes and positive weights required;
%   very large T may underflow and is rejected, not silently altered.
%   The trace phase moments sum(z_j^p)/N are time-independent; the
%   cyclic e1 moments sum(w_j(t)*z_j^p) vary. No continuum density is
%   asserted for this finite unitary model.
%
%   [159] I. Nenciu, IMRN 2005, 647--686, Theorem 3:
%   https://doi.org/10.1155/IMRN.2005.647
%   [160] M. Gekhtman and I. Nenciu, Proc. AMS 136 (2008), 2785--2792:
%   https://arxiv.org/abs/0706.2428
%   [5] B. Simon, OPUC: https://doi.org/10.1090/coll054.1
%   [5] OPUC part 2: https://doi.org/10.1090/coll/054.2
properties={'square','unitary','normal','scalable', ...
    'parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','finite','>=',4,'<=',128});
if mod(n,2)~=0
    error('abl_ladik_cmv:Order','N must be even in this finite boundary model.');
end
if nargin<2||isempty(t),t=0;end
validateattributes(t,{'numeric'},{'scalar','real','finite'});
if nargin<3||isempty(nodes)
    nodes=exp(2i*pi*((0:n-1)'+.5)/n);
end
if nargin<4||isempty(weights),weights=(1:n)';end
if nargin<5||isempty(harmonic),harmonic=1;end
if nargin<6||isempty(output),output='matrix';end
validateattributes(harmonic,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',8});
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('abl_ladik_cmv:Mode','OUTPUT must be a scalar text value.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','operator','parameters','spectral','lax'}))
    error('abl_ladik_cmv:Mode','Invalid output mode.');
end
% Validate nodes, distinctness, and weights using the existing
% inverse_spectral API WITHOUT requiring a potentially ill-conditioned
% coefficient reconstruction until the final requested matrix mode.
initial=inverse_opuc(n,nodes,weights,'cmv','spectral');
z=initial.finite.eigenvalues();
w0=weights(:);w0=w0/sum(w0);
beta=(-1)^(n+1)*conj(prod(z));
if abs(beta+1)>1e-9
    error('abl_ladik_cmv:Boundary', ...
        'Input eigenphases must imply terminal Verblunsky beta=-1.');
end
g=-2*real(z.^harmonic);
logw=log(w0)+t*g;
logw=logw-max(logw);
w=exp(logw);w=w/sum(w);
if any(~isfinite(w))||any(w<=0)
    error('abl_ladik_cmv:WeightUnderflow', ...
        'Positive cyclic weights have underflowed; reduce |T|.');
end
if strcmp(output,'spectral')
    U=inverse_opuc(n,z,w,'cmv','spectral');
    U.family='abl_ladik_cmv';
    U.flow=struct('time',t,'harmonic',harmonic, ...
        'initial_weights',w0,'weights',w, ...
        'generator_values',g,'terminal_alpha',-1, ...
        'hamiltonian','2 Im Tr(C^k)/k');
    U.reference=['Nenciu, IMRN (2005), Thm 3, Eq 4.7; ' ...
        'Gekhtman--Nenciu, Proc AMS (2008), Eq 5.7-5.8'];
    U.notes=[U.notes ' Exact finite AL isospectral Lax solution, ' ...
        'with terminal Verblunsky -1 and fixed determinant. ' ...
        'The cyclic measure varies but finite TRACE weights are 1/N.'];
    return
end
if strcmp(output,'operator')
    U=inverse_opuc(n,z,w,'cmv','operator');return
end
par=inverse_opuc(n,z,w,'cmv','parameters');
if strcmp(output,'parameters')
    U=par;
    U.family='abl_ladik_cmv';
    U.time=t;
    U.harmonic=harmonic;
    U.initial_weights=w0;
    U.generator_values=g;
    U.weight_derivative=w.*(g-sum(w.*g));
    return
end
C=full(inverse_opuc(n,z,w,'cmv'));
if strcmp(output,'matrix'),U=C;return;end
P=C^harmonic;
Pplus=triu(P,1)+diag(diag(P)/2);
B=Pplus-Pplus';
U=struct('matrix',C,'B',B,'commutator',C*B-B*C, ...
    'hamiltonian','2 Im Tr(C^k)/k','time',t, ...
    'harmonic',harmonic,'nodes',z,'cyclic_weights',w, ...
    'initial_weights',w0, ...
    'weight_derivative',w.*(g-sum(w.*g)), ...
    'parameters',par);
end
