function [U,properties,A,B] = unitary_delay_pencil(n,delays,S,mode)
%UNITARY_DELAY_PENCIL Unitary lift of a lossless multichannel delay NEP.
%   U=UNITARY_DELAY_PENCIL(N,DELAYS,S) returns a unitary matrix of order
%   sum(DELAYS), for positive integer delay lengths DELAYS (N entries)
%   and an N-by-N unitary scattering matrix S.
%   The nonlinear secular equation is
%       det(I-diag(z.^DELAYS)*S)=0.
%   These are exactly the eigenvalues of U (including multiplicity).
%   [U,PROPERTIES,A,B] returns a strict pencil A-z*B with
%   B\A=U and A*A'=B*B'. A,B are independently unitary.
%   Defaults: DELAYS=(1:N)', S=DISCRETE_FOURIER(N).
%   The integer-delay finite lift is a benchmark realization of
%   lossless quantum-graph bond scattering/secular equations.
%   [58] https://doi.org/10.1006/aphy.1999.5904
properties={'square','unitary','normal','sparse','scalable','parameter-dependent'};
if nargin==0,U=[];A=[];B=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(delays),delays=(1:n)';end
validateattributes(delays,{'numeric'},{'vector','integer','positive','finite'});
if numel(delays)~=n
    error('unitary_delay_pencil:Delays','DELAYS must have N entries.');
end
delays=double(delays(:));
N=sum(delays);
if ~isfinite(N) || N>flintmax
    error('unitary_delay_pencil:Order','Lift dimension must be finite and exact.');
end
if nargin<3 || isempty(S),S=[];end
if ~isempty(S)
    validateattributes(S,{'numeric'},{'2d','finite','size',[n n]});
    if norm(full(S'*S-eye(n)),'fro')>1e-10*max(1,n)
        error('unitary_delay_pencil:Scattering','S must be unitary.');
    end
end
if nargin<4 || isempty(mode),mode='matrix';end
if strcmpi(char(mode),'spectral')
    U=ep_spectral_reference('unitary_delay_pencil',n,delays);
    A=[];B=[];return
end
if ~strcmpi(char(mode),'matrix')
    error('unitary_delay_pencil:Mode','Supported modes: matrix, spectral.');
end
if isempty(S),S=discrete_fourier(n);end
cuts=[0;cumsum(delays)];
terminal=cuts(2:end);
rows=1:N;
for j=1:n
    idx=cuts(j)+(1:delays(j));
    rows(idx)=idx([2:end 1]);
end
P=sparse(rows,1:N,1,N,N);
C=speye(N);
C(terminal,terminal)=S';
A=C;B=P';
U=P*C;
end
