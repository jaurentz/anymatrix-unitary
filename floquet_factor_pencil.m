function [U,properties,A,B] = floquet_factor_pencil(n,theta1,theta2,mode)
%FLOQUET_FACTOR_PENCIL Strict orthogonal pencil for a split-step quantum walk.
%   The periodic two-state walk is U=Sminus*C2*Splus*C1;
%   A=Splus*C1 and B=(Sminus*C2)' yield B\A=U.
%   The default parameters and 2N-state ordering are exactly those of
%   SPLIT_STEP_WALK. Both pencil coefficients are sparse orthogonal.
%   Kitagawa et al., Phys. Rev. A 82 (2010), 033429.
%   [16] https://doi.org/10.1103/PhysRevA.82.033429
properties={'square','unitary','normal','orthogonal','real','sparse', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];A=[];B=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(theta1),theta1=.3;end
if nargin<3 || isempty(theta2),theta2=.7;end
validateattributes(theta1,{'numeric'},{'scalar','real','finite'});
validateattributes(theta2,{'numeric'},{'scalar','real','finite'});
if nargin<4 || isempty(mode),mode='matrix';end
if strcmpi(char(mode),'spectral')
    U=ep_spectral_reference('floquet_factor_pencil',n,theta1,theta2);
    A=[];B=[];return
end
if ~strcmpi(char(mode),'matrix')
    error('floquet_factor_pencil:Mode','Supported modes: matrix, spectral.');
end
cols=1:2*n;rowplus=cols;rowminus=cols;
rowplus(1:2:end)=2*mod((0:n-1)+1,n)+1;
rowminus(2:2:end)=2*mod((0:n-1)-1,n)+2;
Splus=sparse(rowplus,cols,1,2*n,2*n);
Sminus=sparse(rowminus,cols,1,2*n,2*n);
R1=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
R2=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
C1=kron(speye(n),sparse(R1));
C2=kron(speye(n),sparse(R2));
A=Splus*C1;F=Sminus*C2;B=F';
U=F*A;
end
