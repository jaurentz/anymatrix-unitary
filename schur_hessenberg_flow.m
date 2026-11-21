function [U,properties]=schur_hessenberg_flow(n,t,nodes,weights,representation,method,output)
%SCHUR_HESSENBERG_FLOW QR-factorized isospectral unitary Hessenberg flow.
%   U=SCHUR_HESSENBERG_FLOW(N,T,NODES,WEIGHTS,REP,METHOD,OUTPUT)
%   starts with the upper-Hessenberg realization H0 constructed from
%   distinct circle nodes z_j and positive e1 weights w_j(0). For
%   real T its positive-diagonal QR factorization defines
%
%       exp(T*H0) = Q(T)*R(T),  diag(R(T)) > 0,
%       H(T) = Q(T)'*H0*Q(T).
%
%   This is unitary upper-Hessenberg, isospectral to H0 in exact
%   arithmetic. Its e1 weights become
%
%       w_j(T) = w_j(0)*exp(2*T*real(z_j)) / Z(T).
%
%   REP='hessenberg' (default) uses the direct QR path (METHOD='qr')
%   or an independent inverse-OPUC reconstruction (METHOD='inverse').
%   REP='cmv' uses the same evolving nodes and cyclic weights via
%   inverse_opuc(...,'cmv'); it is a CROSS-REPRESENTATION of the flow,
%   not itself the Hessenberg QR similarity.
%
%   OUTPUT='matrix' (default), 'parameters', 'operator', or 'spectral'.
%   'operator' deliberately uses a core-product inverse reconstruction,
%   avoiding the dense exponential and O(N^3) QR of 'matrix'.
%   The finite trace spectrum/moments are constant; the e1 cyclic
%   spectral weights evolve. Matrix generation can lose accuracy for
%   large |T| when exponential QR becomes rank-deficient; spectral mode
%   uses centered log weights and avoids that factorization entirely.
%
%   Ammar--Gragg, Schur Flows for Orthogonal Hessenberg Matrices (1993):
%   https://calhoun.nps.edu/bitstream/handle/10945/30447/Ammar_Gragg_Schur.pdf?sequence=3
%   Golinskii, Sb. Math. 197 (2006):
%   https://doi.org/10.1070/SM2006v197n08ABEH003792
%   Simon, OPUC (2005): https://doi.org/10.1090/coll054.1
%   Simon, OPUC part 2: https://doi.org/10.1090/coll/054.2
properties={'square','unitary','normal','scalable', ...
    'parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite','<=',512});
if nargin<2||isempty(t),t=0;end
validateattributes(t,{'numeric'},{'scalar','real','finite'});
if nargin<3||isempty(nodes)
    nodes=exp(2i*pi*((0:n-1)'+.17)/n);
end
if nargin<4||isempty(weights),weights=ones(n,1)/n;end
if nargin<5||isempty(representation),representation='hessenberg';end
if nargin<6||isempty(method),method='qr';end
if nargin<7||isempty(output),output='matrix';end
if ~(ischar(representation)||(isstring(representation)&&isscalar(representation))) ...
        ||~(ischar(method)||(isstring(method)&&isscalar(method))) ...
        ||~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('schur_hessenberg_flow:Options','REP, METHOD, OUTPUT must be strings.');
end
representation=lower(char(representation));
method=lower(char(method));
output=lower(char(output));
if ~any(strcmp(representation,{'hessenberg','cmv'})) ...
        ||~any(strcmp(method,{'qr','inverse'})) ...
        ||~any(strcmp(output,{'matrix','operator','parameters','spectral'}))
    error('schur_hessenberg_flow:Options','Invalid representation, method or output.');
end
initial=inverse_opuc(n,nodes,weights,'hessenberg','spectral');
z=initial.finite.eigenvalues();
w0=weights(:);w0=w0/sum(w0);
logw=log(w0)+2*t*real(z);
logw=logw-max(logw);
w=exp(logw);
w=w/sum(w);
if any(w==0)||any(~isfinite(w))
    error('schur_hessenberg_flow:WeightUnderflow', ...
        'Cyclic weights underflowed; reduce |T|.');
end
if strcmp(output,'spectral')
    U=inverse_opuc(n,z,w,representation,'spectral');
    U.family='schur_hessenberg_flow';
    U.flow=struct('time',t,'nodes',z, ...
        'initial_weights',w0,'weights',w, ...
        'representation',representation,'method',method);
    U.reference='Ammar--Gragg (1993); Golinskii (2006); Simon (2005)';
    U.notes=[U.notes ' Positive-diagonal Hessenberg QR flow keeps ' ...
        'equal-weight trace phases fixed and changes e1 weights. ' ...
        'CMV is an inverse-spectral cross-representation.'];
    return
end
if strcmp(output,'operator')
    U=inverse_opuc(n,z,w,representation,'operator');
    return
end
if strcmp(representation,'cmv')||strcmp(method,'inverse')
    if strcmp(output,'parameters')
        U=inverse_opuc(n,z,w,representation,'parameters');
        U.time=t;U.initial_weights=w0;
        U.method='inverse_spectral';
    else
        U=inverse_opuc(n,z,w,representation,'matrix');
    end
    return
end
H0=full(inverse_opuc(n,z,w0,'hessenberg','matrix'));
% Scaling by a positive scalar does not alter QR's positive-diagonal Q.
% Re(z) <= 1, so exp(t*H0-|t|I) does not overflow spectrally.
E=expm(t*H0-abs(t)*eye(n));
[Q,R]=qr(E);
d=diag(R);
if any(abs(d)<=eps(max(1,norm(E,'fro')))*n)
    error('schur_hessenberg_flow:QRBreakdown', ...
        'QR exponential is numerically rank deficient; use inverse mode.');
end
Q=Q*diag(d./abs(d));  % positive-diagonal QR convention
H=Q'*H0*Q;
if strcmp(output,'parameters')
    U=struct('initial_hessenberg',H0,'hessenberg',H, ...
        'Q',Q,'time',t,'nodes',z, ...
        'initial_weights',w0,'weights',w, ...
        'representation',representation,'method','qr', ...
        'hessenberg_residual',norm(tril(H,-2),'fro'));
else
    U=H;
end
end
