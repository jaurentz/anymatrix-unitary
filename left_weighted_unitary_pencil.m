function [U,properties] = left_weighted_unitary_pencil(n,kappa,mode)
%LEFT_WEIGHTED_UNITARY_PENCIL Strict unitary pencil with ill-scaled factors.
%   For P=CYCLE_SHIFT(N), B=diag(logspace(0,log10(KAPPA),N)),
%   A=B*P. Thus B\A=P is exactly unitary and AA'=BB', even though
%   A and B are not themselves unitary when KAPPA~=1.
%   Use to test generalized eigensolvers and coefficient scaling.
%   KAPPA>=1 (default 1e3). This is a synthetic test construction,
%   based on structured eigenproblem principles, not an attributed model.
%   [4] https://doi.org/10.1137/1.9781611975345
%   Terminal 'pencil' mode returns struct P with P.A, P.B, P.U.
%   Standard [U,PROPERTIES] signature preserves Anymatrix discovery.
properties={'square','unitary','orthogonal','real','normal','sparse', ...
    'scalable','eigenvalues','parameter-dependent'};
if nargin==0,U=[];A=[];B=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(kappa),kappa=1e3;end
validateattributes(kappa,{'numeric'},{'scalar','real','finite','>=',1});
if nargin<3 || isempty(mode),mode='matrix';end
if strcmpi(char(mode),'spectral')
    U=ep_spectral_reference('left_weighted_unitary_pencil',n);
    A=[];B=[];return
end
if ~any(strcmpi(char(mode),{'matrix','pencil'}))
    error('left_weighted_unitary_pencil:Mode','Supported modes: matrix, pencil, spectral.');
end
d=exp(linspace(0,log(kappa),n).');
B=spdiags(d,0,n,n);
U=cycle_shift(n);
A=B*U;
if strcmpi(char(mode),'pencil')
    U=struct('A',A,'B',B,'U',U);
end
end
