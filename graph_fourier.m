function [U,properties]=graph_fourier(n,adjacency,output)
%GRAPH_FOURIER Orthogonal graph Fourier analysis from a weighted Laplacian.
%   U=GRAPH_FOURIER(N,A) uses an N-by-N undirected symmetric nonnegative
%   adjacency matrix A with zero diagonal. Default A is the path graph.
%   L=diag(sum(A,2))-A=V*diag(lambda)*V'; U=V'.
%   Columns of V are in ascending Laplacian-frequency order.
%   Within repeated eigenspaces, projected standard basis vectors are
%   orthonormalized in vertex order to fix a deterministic convention.
%   The eigenvalues of L are NOT the eigenvalues of U.
%   Only an explicit matrix is returned; a dense cached eigensystem
%   is not a matrix-free implementation.
%   Reference: D. I. Shuman et al., IEEE SPM 30 (2013).
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [57] https://doi.org/10.1109/MSP.2012.2235192
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('graph_fourier',n,adjacency);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2 || isempty(adjacency)
    if n==1
        A=sparse(1,1);
    else
        A=sparse(1:n-1,2:n,1,n,n);
        A=A+A';
    end
else
    A=adjacency;
end
validateattributes(A,{'numeric'}, ...
    {'size',[n n],'real','finite','nonnegative'});
if norm(A-A','fro')>1e-12*max(1,norm(A,'fro')) || any(diag(A)~=0)
    error('graph_fourier:Adjacency', ...
        'Adjacency must be symmetric with zero diagonal.');
end
L=diag(full(sum(A,2)))-full(A);
[V,D]=eig(L);
[d,idx]=sort(real(diag(D)),'ascend');
V=real(V(:,idx));
tol=1e-10*max(1,max(abs(d)));
first=1;
while first<=n
    last=first;
    while last<n && abs(d(last+1)-d(first))<=tol
        last=last+1;
    end
    block=V(:,first:last);
    count=last-first+1;
    if count>1
        P=block*block';
        Q=zeros(n,count);
        used=0;
        for j=1:n
            v=P(:,j);
            for pass=1:2
                if used>0,v=v-Q(:,1:used)*(Q(:,1:used)'*v);end
            end
            r=norm(v);
            if r>1e-9
                used=used+1;
                Q(:,used)=v/r;
                if used==count,break;end
            end
        end
        if used~=count,error('graph_fourier:Basis','Degenerate subspace canonicalization failed.');end
        V(:,first:last)=Q;
    end
    for j=first:last
        [~,pivot]=max(abs(V(:,j)));
        if V(pivot,j)<0,V(:,j)=-V(:,j);end
    end
    first=last+1;
end
U=V';
end
