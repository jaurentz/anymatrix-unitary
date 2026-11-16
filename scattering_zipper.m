function [Z,properties] = scattering_zipper(n,q,blocks,leftBoundary,rightBoundary,mode)
%SCATTERING_ZIPPER General finite multichannel unitary scattering zipper.
%   Z=SCATTERING_ZIPPER(N,Q) gives an N-site, Q-channel sparse unitary
%   of order N*Q, N even >=2. Default Q=2, with nontrivial channel mixers.
%   BLOCKS is 2Q-by-2Q-by-(N-1) with S_2,...,S_N; each S_j must be
%   unitary with invertible upper-right Q-by-Q block. LEFTBOUNDARY and
%   RIGHTBOUNDARY are independent Q-by-Q unitaries (defaults I).
%   Z=V*W, V=diag(S_2,S_4,...,S_N),
%   W=diag(LEFTBOUNDARY,S_3,S_5,...,S_(N-1),RIGHTBOUNDARY).
%   'pencil' returns P.A=W, P.B=V', P.U=Z, with P.B\P.A=Z.
%   'spectral' returns an unknown-trace-spectrum descriptor.
%   Arbitrary U(2Q) scatterers are more general than block-CMV cores.
%   Marin and Schulz-Baldes, J. Spectr. Theory 3 (2013), 47-82.
%   https://doi.org/10.4171/JST/37
properties={'square','unitary','normal','sparse','scalable','parameter-dependent'};
if nargin==0,Z=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','finite','>=',2});
if mod(n,2),error('scattering_zipper:EvenSites','N must be even.');end
if nargin<2 || isempty(q),q=2;end
validateattributes(q,{'numeric'},{'scalar','integer','finite','positive'});
if nargin<3,blocks=[];end
if nargin<4 || isempty(leftBoundary),leftBoundary=eye(q);end
if nargin<5 || isempty(rightBoundary),rightBoundary=eye(q);end
if nargin<6 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('scattering_zipper:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if strcmp(mode,'spectral')
    Z=unitary_spectral_reference('finite_cmv',n*q);
    Z.family='scattering_zipper';
    Z.reference='Marin and Schulz-Baldes, J. Spectr. Theory 3 (2013)';
    Z.notes=['General multichannel zipper with unitary endpoint scattering; ' ...
        'no universal finite or limiting trace spectrum is asserted.'];
    return
end
if ~any(strcmp(mode,{'matrix','pencil'}))
    error('scattering_zipper:Mode','Supported modes: matrix, pencil, spectral.');
end
if isempty(blocks)
    F=fft(eye(q))/sqrt(q);
    blocks=complex(zeros(2*q,2*q,n-1));
    for j=2:n
        a=diag(.25*exp(1i*(.17*j+.23*(1:q))));
        rho=sqrtm(eye(q)-a*a'); rt=sqrtm(eye(q)-a'*a);
        u=diag(exp(1i*.13*j*(1:q)))*F;
        v=F'*diag(exp(-1i*.19*j*(1:q)));
        blocks(:,:,j-1)=[a,rho*u;v*rt,-v*a'*u];
    end
end
if ~(isnumeric(blocks)&&all(isfinite(blocks(:))) ...
    &&size(blocks,1)==2*q&&size(blocks,2)==2*q ...
    &&size(blocks,3)==n-1)
    error('scattering_zipper:Blocks','BLOCKS must be 2Q-by-2Q-by-(N-1).');
end
validateattributes(leftBoundary,{'numeric'},{'2d','finite','size',[q q]});
validateattributes(rightBoundary,{'numeric'},{'2d','finite','size',[q q]});
tol=1e-10*max(q,1);
if norm(leftBoundary'*leftBoundary-eye(q),'fro')>tol ...
    ||norm(rightBoundary'*rightBoundary-eye(q),'fro')>tol
    error('scattering_zipper:Boundary','Endpoint scatterers must be unitary.');
end
for j=1:n-1
    s=blocks(:,:,j);
    if norm(s'*s-eye(2*q),'fro')>1e-10*max(1,2*q) ...
        ||rcond(s(1:q,q+1:2*q))<=eps
        error('scattering_zipper:Blocks', ...
            'Scatterers must be unitary with invertible upper-right block.');
    end
end
N=n*q;
V=speye(N);W=speye(N);
W(1:q,1:q)=leftBoundary;
W(N-q+1:N,N-q+1:N)=rightBoundary;
for j=2:n
    idx=(j-2)*q+(1:2*q);
    if mod(j,2)==0,V(idx,idx)=blocks(:,:,j-1);
    else,W(idx,idx)=blocks(:,:,j-1);end
end
Z=V*W;
if strcmp(mode,'pencil')
    Z=struct('A',W,'B',V','U',Z);
end
end
