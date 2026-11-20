function [U,properties]=optical_mesh(n,layout,theta,phi,external,output)
%OPTICAL_MESH Lossless Reck triangular or Clements rectangular mesh.
%   U=OPTICAL_MESH(N,LAYOUT,THETA,PHI,EXTERNAL,OUTPUT) has N
%   optical modes and N*(N-1)/2 ordered adjacent-mode couplers.
%   LAYOUT='reck' uses a triangular sequential elimination pattern;
%   LAYOUT='clements' uses N alternating disjoint nearest-neighbor
%   layers. One 2-by-2 core acts on modes j,j+1 as
%      G=[exp(1i*phi)*cos(theta), -sin(theta);
%         exp(1i*phi)*sin(theta),  cos(theta)],
%   with 0<=theta<=pi/2 and real phi. Couplers are applied in the
%   order returned by the chosen geometry: U=D*G_M*...*G_1,
%   D=diag(exp(1i*EXTERNAL)). THETA and PHI must each have M entries.
%   Default parameters give deterministic nontrivial unitary meshes.
%   OUTPUT='operator'/'implicit' returns an O(N^2) local-core action
%   with adjoint and multiple-RHS support; 'spectral' gives an exact
%   finite spectrum only in the diagonal (THETA=0) case and otherwise
%   labels the general analytical trace eigenphase law unknown.
%   No optical losses or nonunitary contractions are modeled.
%   [118] https://doi.org/10.1103/PhysRevLett.73.58
%   [119] https://doi.org/10.1364/OPTICA.3.001460
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(layout),layout='clements';end
if ~(ischar(layout)||(isstring(layout)&&isscalar(layout)))
    error('optical_mesh:Layout','LAYOUT must be reck or clements.');
end
layout=lower(char(layout));
switch layout
    case 'reck'
        pairs=zeros(n*(n-1)/2,1);t=0;
        for col=1:n-1
            for row=n-1:-1:col
                t=t+1;pairs(t)=row;
            end
        end
    case 'clements'
        pairs=zeros(n*(n-1)/2,1);t=0;
        for layer=1:n
            start=1+mod(layer-1,2);
            for row=start:2:n-1
                t=t+1;pairs(t)=row;
            end
        end
    otherwise
        error('optical_mesh:Layout','LAYOUT must be reck or clements.');
end
m=numel(pairs);
if nargin<3||isempty(theta),theta=pi/6+.11*sin((1:m)'*.7);end
if nargin<4||isempty(phi),phi=.25*(1:m)';end
if nargin<5||isempty(external),external=.17*(0:n-1)';end
if ~(isnumeric(theta)&&isreal(theta)&&isvector(theta) ...
        &&numel(theta)==m&&all(isfinite(theta(:))) ...
        &&all(theta(:)>=0)&&all(theta(:)<=pi/2))
    error('optical_mesh:Theta','THETA must contain M angles in [0,pi/2].');
end
if ~(isnumeric(phi)&&isreal(phi)&&isvector(phi) ...
        &&numel(phi)==m&&all(isfinite(phi(:))))
    error('optical_mesh:Phi','PHI must contain M finite real phases.');
end
if ~(isnumeric(external)&&isreal(external)&&isvector(external) ...
        &&numel(external)==n&&all(isfinite(external(:))))
    error('optical_mesh:External','EXTERNAL must have N real finite phases.');
end
theta=theta(:);phi=phi(:);external=external(:);
if nargin<6||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('optical_mesh:Mode','Invalid output mode.');
end
c=cos(theta);s=sin(theta);phase=exp(1i*phi);
ext=exp(1i*external);
F=@(X,flag)um_optical_mesh_apply(X,flag,n,pairs,c,s,phase,ext);
switch lower(char(output))
    case {'operator','implicit'}
        U=F;
    case {'matrix','explicit'}
        U=F(eye(n),'notransp');
    case 'spectral'
        U=unitary_spectral_reference('optical_mesh',n);
        U.reference='Reck et al. PRL 73 (1994); Clements et al. Optica 3 (2016)';
        if all(theta==0)
            z=ext;
            for j=1:m
                z(pairs(j))=z(pairs(j))*phase(j);
            end
            U.finite.kind='exact';
            U.finite.eigenvalues=@()z;
            U.finite.moments=@(k)um_optical_mesh_moments(k,z);
            U.notes='Diagonal zero-coupling network: exact finite eigenphase multiset.';
        else
            U.notes=['General product of lossless photonic beamsplitters: ' ...
                'no claimed closed-form trace eigenphase law.'];
        end
    otherwise
        error('optical_mesh:Mode','Use matrix, operator, or spectral.');
end
end
