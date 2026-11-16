function [U,properties,A,B] = finite_cmv_pencil(n,alpha,mode)
%FINITE_CMV_PENCIL Strict unitary pencil of a scalar finite CMV matrix.
%   U=FINITE_CMV_PENCIL(N,ALPHA) equals FINITE_CMV(N,ALPHA).
%   ALPHA is length N, with interior moduli <1 and terminal modulus 1.
%   [U,PROPERTIES,A,B] returns A=M, B=L' from the CMV factorization
%   U=L*M, so B\A=U. Both A and B are sparse unitary matrices.
%   [6] https://doi.org/10.1016/S0024-3795(02)00457-3
%   [7] https://doi.org/10.1016/j.cam.2006.10.033
properties={'square','unitary','normal','sparse','banded','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];A=[];B=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(alpha),alpha=[zeros(n-1,1);1];end
validateattributes(alpha,{'numeric'},{'vector','finite'});
if numel(alpha)~=n
    error('finite_cmv_pencil:Alpha','ALPHA must have N entries.');
end
alpha=alpha(:);
if any(abs(alpha(1:end-1))>=1) || ...
        abs(abs(alpha(end))-1)>128*eps
    error('finite_cmv_pencil:Alpha','Interior coefficients must be contractive and closure unimodular.');
end
if nargin<3 || isempty(mode),mode='matrix';end
if strcmpi(char(mode),'spectral')
    U=ep_spectral_reference('finite_cmv_pencil',n);
    A=[];B=[];return
end
if ~strcmpi(char(mode),'matrix')
    error('finite_cmv_pencil:Mode','Supported modes: matrix, spectral.');
end
[U,~,A,B]=block_cmv_pencil(n,1,reshape(alpha(1:end-1),1,1,[]),alpha(end));
end
