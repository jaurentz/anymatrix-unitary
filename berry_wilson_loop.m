function [U,properties]=berry_wilson_loop(n,frames,output)
%BERRY_WILSON_LOOP Discrete Berry/Wilczek-Zee holonomy from moving frames.
%   W=BERRY_WILSON_LOOP(N,FRAMES) uses M>=2 orthonormal D-by-N frames
%   FRAMES(:,:,j) on a CLOSED parameter loop. The final edge connects
%   frame M back to frame 1; DO NOT repeat the initial frame as an
%   additional last sample. N is the occupied rank and matrix order.
%   The raw overlaps M_j=Q_j'*Q_(j+1) are usually contractions.
%   Provided all overlaps are nonsingular, their UNITARY polar factors
%      V_j=polar(M_j)=U_j*Z_j' for svd(M_j)=U_j*S_j*Z_j'
%   give the strictly unitary W=V_1*...*V_M.
%   Under independent frame changes Q_j -> Q_j*G_j, each link
%   transforms V_j -> G_j'*V_j*G_(j+1) and W -> G_1'*W*G_1;
%   only the eigenphases and conjugacy invariants are gauge invariant.
%
%   The default FRAMES are a deterministic 2-state latitude loop for
%   N=1, or a closed D=N+2, noncommuting rank-N loop for N>1.
%   'operator' uses existing UM_WILSON_APPLY with O(M*N^2) actions,
%   reverse-order adjoint and complex multiple RHS, no full W matrix.
%   'links' exposes the computed unitary polar overlap links.
%   'data' returns raw overlaps, singular values, links, their raw
%   product and exact sampled holonomy, to distinguish contractions
%   from strict finite unitary matrices.
%   'spectral' is unknown for generic non-Abelian loops. For N=1
%   it gives the EXACT sampled polar Berry phase and integer moments,
%   not an asserted topological bulk invariant or continuum limit.
%   Singular or nearly singular overlaps are rejected; numerical polar
%   factors need a well-conditioned subspace overlap.
%   For optical/Bloch band data, input frames must actually span the
%   intended physical occupied eigenspaces. No hidden band solver is used.
%
%   [124] Wilczek and Zee, PRL 52 (1984), 2111:
%   https://doi.org/10.1103/PhysRevLett.52.2111
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(frames),frames=default_frames(n);end
if ~isnumeric(frames)||ndims(frames)>3||size(frames,2)~=n ...
        ||size(frames,1)<n||size(frames,3)<2 ...
        ||any(~isfinite(frames(:)))
    error('berry_wilson_loop:Frames', ...
        'FRAMES must be finite D-by-N-by-M (D>=N,M>=2).');
end
frames=double(frames);
d=size(frames,1);
m=size(frames,3);
for j=1:m
    Q=frames(:,:,j);
    if norm(Q'*Q-eye(n),'fro')>2e-10*max(1,n)
        error('berry_wilson_loop:Orthonormal', ...
            'Each occupied-band frame must have orthonormal columns.');
    end
end
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('berry_wilson_loop:Mode','OUTPUT must be a scalar string.');
end
mode=lower(char(output));
if ~any(strcmp(mode,{'matrix','explicit','operator','implicit', ...
        'spectral','links','data'}))
    error('berry_wilson_loop:Mode', ...
        'Use matrix, operator, spectral, links, or data.');
end
overlaps=complex(zeros(n,n,m));
links=complex(zeros(n,n,m));
singular_values=zeros(n,m);
for j=1:m
    next=mod(j,m)+1;
    overlap=frames(:,:,j)'*frames(:,:,next);
    [P,S,Q]=svd(overlap);
    sigma=diag(S);
    if min(sigma)<=1e-8
        error('berry_wilson_loop:SingularOverlap', ...
            'A frame overlap is singular or ill conditioned.');
    end
    overlaps(:,:,j)=overlap;
    links(:,:,j)=P*Q';
    singular_values(:,j)=sigma;
end
switch mode
    case {'matrix','explicit'}
        U=lattice_wilson_loop(n,links,'matrix');
    case {'operator','implicit'}
        U=lattice_wilson_loop(n,links,'operator');
    case 'links'
        U=links;
    case 'data'
        raw=eye(n);
        for j=1:m,raw=raw*overlaps(:,:,j);end
        U=struct('dimension',n,'ambient_dimension',d,'samples',m, ...
            'overlaps',overlaps,'singular_values',singular_values, ...
            'links',links,'raw_product',raw, ...
            'holonomy',lattice_wilson_loop(n,links,'matrix'), ...
            'notes','Numerical finite polar-overlap holonomy, not a bulk invariant.');
    case 'spectral'
        U=unitary_spectral_reference('berry_wilson_loop',n);
        U.reference='Wilczek and Zee (1984), doi:10.1103/PhysRevLett.52.2111';
        if n==1
            phase=prod(reshape(links(1,1,:),[],1));
            U.finite.kind='exact';
            U.finite.eigenvalues=@()phase;
            U.finite.moments=@(k)local_moments(phase,k);
            U.finite.atoms=struct('angles',angle(phase),'weights',1);
            U.notes=['Exact finite Abelian polar-overlap phase of the ' ...
                'supplied discrete loop; no inferred continuum winding number.'];
        else
            U.notes=['Finite non-Abelian holonomy computed by polar ' ...
                'overlaps; generic analytical trace eigenphases unknown.'];
        end
end
end

function frames=default_frames(n)
% Closed and reproducible. N=1 has a latitude spinor with known
% continuum Berry phase; N>1 uses noncommuting Grassmannian frames.
if n==1
    theta=1.1;m=32;
    frames=complex(zeros(2,1,m));
    for j=1:m
        t=2*pi*(j-1)/m;
        frames(:,1,j)=[cos(theta/2);exp(1i*t)*sin(theta/2)];
    end
else
    d=n+2;m=16;
    Q0=fft(eye(d))/sqrt(d);
    Q0=Q0(:,1:n);
    B=complex(zeros(d));
    for j=1:d-1
        B(j,j+1)=.38+.07i*j;
        B(j+1,j)=-conj(B(j,j+1));
    end
    frames=complex(zeros(d,n,m));
    for j=1:m
        t=2*pi*(j-1)/m;
        D=diag(exp(1i*(0:d-1)*t));
        frames(:,:,j)=D*expm(sin(t)*B)*Q0;
    end
end
end

function out=local_moments(z,k)
if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
        ||any(k(:)~=round(k(:)))
    error('berry_wilson_loop:Moments','Moment orders must be integers.');
end
out=z.^k;
end
