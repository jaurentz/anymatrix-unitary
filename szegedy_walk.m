function [U,properties]=szegedy_walk(n,P,mode)
%SZEGEDY_WALK Product-of-reflections quantum walk for a Markov chain.
%   U=SZEGEDY_WALK(N,P) constructs an N^2-by-N^2 unitary matrix.
%   P is an N-by-N row-stochastic matrix with nonnegative entries.
%   State |i,j> has index (i-1)*N+j. A and B are N^2-by-N isometries:
%      A(:,i)=sum_j sqrt(P(i,j)) |i,j>,
%      B(:,j)=sum_i sqrt(P(j,i)) |i,j> = SWAP*A(:,j).
%   U=(2*B*B'-I)*(2*A*A'-I).
%   With D=A'*B=sqrt(P.*P.'), the nontrivial eigenvalues follow
%   exp(+/-2i*acos(svd(D))) with correct 0/1 and complement phases.
%   N^2 may be large: MODE='operator' provides an O(N^2) matvec
%   without assembling the generally dense matrix; 'spectral' returns
%   a lazy exact finite spectral descriptor using just the N-by-N D.
%   Default P is a lazy nearest-neighbor Markov chain on a cycle.
%
%   Mario Szegedy, Quantum speed-up of Markov chain based algorithms,
%   FOCS 2004, 32-41. https://doi.org/10.1109/FOCS.2004.53
properties={'square','unitary','normal','real','scalable','eigenvalues', ...
            'parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(P)
    if n==1
        P=1;
    else
        P=sparse(n,n);
        for i=1:n
            P(i,i)=.5;
            j=mod(i,n)+1; h=mod(i-2,n)+1;
            P(i,j)=P(i,j)+.25;
            P(i,h)=P(i,h)+.25;
        end
    end
end
validateattributes(P,{'numeric'},{'2d','real','nonnegative','finite','size',[n n]});
if max(abs(full(sum(P,2))-1))>1e-12
    error('szegedy_walk:Stochastic','Each row of P must sum to one.');
end
if nargin<3 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('szegedy_walk:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
N=n*n;
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('finite_cmv',N);
    U.family='szegedy_walk';
    U.reference='Szegedy, FOCS (2004), DOI 10.1109/FOCS.2004.53';
    U.notes=['Exact finite eigenvalues from singular values of the ' ...
        'N-by-N Markov discriminant, including all +/-1 multiplicities. ' ...
        'No universal limiting spectral density is asserted.'];
    U.finite.kind='exact';
    U.finite.eigenvalues=@()szegedy_spectrum(P);
    U.finite.moments=@(k)szegedy_moments(P,k);
    return
end
if strcmp(mode,'operator')
    U=@(X,flag)um_szegedy_apply(X,flag,P);
    return
end
if ~strcmp(mode,'matrix')
    error('szegedy_walk:Mode','Supported modes: matrix, operator, spectral.');
end
if n>64
    error('szegedy_walk:ExplicitOrder', ...
          'Explicit N^2-by-N^2 construction is too costly; use operator mode.');
end
rows=(1:N)';
cols=kron((1:n)',ones(n,1));
vals=reshape(sqrt(full(P)).',[],1);
A=sparse(rows,cols,vals,N,n);
reverse=reshape(reshape(1:N,n,n).',[],1);
S=sparse(reverse,1:N,1,N,N);
B=S*A;
RA=2*(A*A')-speye(N);
RB=2*(B*B')-speye(N);
U=RB*RA;
end

function z=szegedy_spectrum(P)
n=size(P,1);
D=sqrt(full(P).*full(P.'));
s=svd(D);
tol=128*eps(max(n,1));
onesig=s>=1-tol;
zerosig=s<=tol;
middle=~(onesig|zerosig);
angles=2*acos(min(1,max(0,s(middle))));
nminus=2*nnz(zerosig);
nplus=n*n-2*nnz(middle)-nminus;
if nplus<0,error('szegedy_walk:Spectrum','Invalid multiplicity.');end
z=[exp(1i*angles);exp(-1i*angles); ...
   -ones(nminus,1);ones(nplus,1)];
end

function m=szegedy_moments(P,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
z=szegedy_spectrum(P);
m=reshape(mean(bsxfun(@power,z(:),reshape(k,1,[])),1),size(k));
end
