function [U,properties]=gauss_recurrence_transform(n,diagonal,offdiagonal,output)
%GAUSS_RECURRENCE_TRANSFORM Orthogonal transform for a Jacobi recurrence.
%   Q=GAUSS_RECURRENCE_TRANSFORM(N,B,A) is the real orthogonal
%   Gauss-quadrature evaluation transform for the symmetric Jacobi matrix
%   J=diag(B)+diag(A,1)+diag(A,-1). B has N finite real entries,
%   and A has N-1 strictly positive real entries.
%   The eigenvalues of J are quadrature NODES, not eigenphases of Q.
%   Q is the transpose of the ordered J eigenvector matrix, with
%   eigenvectors made positive in their first component. Thus
%   Q*J*Q'=diag(sorted eigenvalues of J) and Q'*Q=I.
%   Default B=0 and A=sqrt(k/2), matching Gauss--Hermite.
%   'spectral' reports finite trace eigenphase law as UNKNOWN; it
%   does not infer eigenphases from the quadrature nodes.
%   Dense Golub--Welsch setup O(N^3), multiply O(N^2). No fast
%   matrix-free action is claimed for arbitrary user recurrences.
%   [104] https://doi.org/10.1090/S0025-5718-69-99647-1
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(diagonal),diagonal=zeros(n,1);end
if nargin<3||isempty(offdiagonal)
    offdiagonal=sqrt((1:n-1)'/2);
end
if ~(isnumeric(diagonal)&&isreal(diagonal)&&isvector(diagonal) ...
        &&numel(diagonal)==n&&all(isfinite(diagonal(:))))
    error('gauss_recurrence_transform:Diagonal', ...
        'B must be a finite real N-vector.');
end
if ~(isnumeric(offdiagonal)&&isreal(offdiagonal) ...
        &&numel(offdiagonal)==n-1 ...
        &&(n==1||isvector(offdiagonal)) ...
        &&all(isfinite(offdiagonal(:)))&&all(offdiagonal(:)>0))
    error('gauss_recurrence_transform:Offdiagonal', ...
        'A must have N-1 finite strictly positive real entries.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('gauss_recurrence_transform:Mode','Invalid output mode.');
end
switch lower(char(output))
    case {'matrix','explicit'}
        U=um_gauss_basis(diagonal(:),offdiagonal(:));
    case 'spectral'
        U=unitary_spectral_reference('gauss_recurrence_transform',n);
        U.reference='Golub and Welsch, Math. Comp. 23 (1969)';
        U.notes=['Jacobi eigenspectrum gives real Gaussian quadrature nodes, ' ...
            'not unit-circle eigenphases of the orthogonal transform; ' ...
            'finite trace spectral law unknown in general.'];
    otherwise
        error('gauss_recurrence_transform:Mode', ...
            'Only matrix and spectral output are supported.');
end
end
