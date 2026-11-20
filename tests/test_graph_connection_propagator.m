function tests=test_graph_connection_propagator
tests=functiontests(localfunctions);
end

function test_default_cycle_exact_phase_formula_and_hermiticity(tc)
for n=[3 4 7 11]
    for q=[1 2 3]
        t=.73;
        D=graph_connection_propagator(n,q,[],[],[],t,'laplacian');
        L=full(D.L);
        N=n*q;
        verifyEqual(tc,size(L),[N N]);
        verifyLessThan(tc,norm(L-L','fro'),2e-13*N);
        verifyGreaterThanOrEqual(tc,min(eig(L)),-2e-11);
        U=graph_connection_propagator(n,q,[],[],[],t);
        verifyEqual(tc,size(U),[N N]);
        verifyLessThan(tc,norm(U'*U-eye(N),'fro'),5e-11*N);
        S=graph_connection_propagator(n,q,[],[],[],t,'spectral');
        verifyEqual(tc,S.finite.kind,'exact');
        verifyEqual(tc,S.dimension,N);
        z=S.finite.eigenvalues();
        verifyEqual(tc,numel(z),N);
        eigs=eig(U);
        for j=1:N
            verifyLessThan(tc,min(abs(z-eigs(j))),4e-9);
        end
        for k=[-5 -2 -1 0 1 2 5]
            verifyLessThan(tc,abs(S.finite.moments(k)-trace(U^k)/N), ...
                6e-10);
        end
        V=graph_connection_propagator(n,q,[],[],[],-t);
        verifyLessThan(tc,norm(V-U','fro'),3e-11*N);
    end
end
end

function test_connection_energy_and_nonabelian_frame_covariance(tc)
n=4;q=2;
edges=[1 2;2 3;3 4;4 1;1 3];
weights=[1;.75;1.2;.8;.4];
m=size(edges,1);
R=complex(zeros(q,q,m));
for j=1:m
    a=.18*j;
    Q=[cos(a),-sin(a);sin(a),cos(a)];
    R(:,:,j)=Q*diag(exp(1i*[.07*j,-.12*j]));
end
D=graph_connection_propagator(n,q,edges,weights,R,1.1,'laplacian');
L=full(D.L);
verifyLessThan(tc,norm(L-L','fro'),3e-13);
verifyGreaterThanOrEqual(tc,min(eig(L)),-2e-11);
x=randn(n*q,1)+1i*randn(n*q,1);
energy=0;
for e=1:m
    ii=(edges(e,1)-1)*q+(1:q);
    jj=(edges(e,2)-1)*q+(1:q);
    energy=energy+weights(e)*norm(x(ii)-R(:,:,e)*x(jj))^2;
end
verifyLessThan(tc,abs(real(x'*L*x)-energy),2e-11);
G=complex(zeros(q,q,n));globalFrame=zeros(n*q);
for i=1:n
    G(:,:,i)=diag(exp(1i*[.13*i,-.21*i]));
    ii=(i-1)*q+(1:q);
    globalFrame(ii,ii)=G(:,:,i);
end
Rg=R;
for e=1:m
    i=edges(e,1);j=edges(e,2);
    Rg(:,:,e)=G(:,:,i)*R(:,:,e)*G(:,:,j)';
end
Lnew=graph_connection_propagator(n,q,edges,weights,Rg,1.1,'laplacian');
verifyLessThan(tc,norm(full(Lnew.L)-globalFrame*L*globalFrame','fro'), ...
    3e-11);
U=graph_connection_propagator(n,q,edges,weights,R,1.1);
Ug=graph_connection_propagator(n,q,edges,weights,Rg,1.1);
verifyLessThan(tc,norm(Ug-globalFrame*U*globalFrame','fro'),4e-10);
S=graph_connection_propagator(n,q,edges,weights,R,1.1,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
end

function test_sparse_cayley_is_unitary_but_not_exact_exponential(tc)
n=5;q=2;t=1.36;
U=graph_connection_propagator(n,q,[],[],[],t);
errors=zeros(1,3);
for j=1:3
    s=4*2^(j-1);
    C=graph_connection_propagator(n,q,[],[],[],t,'cayley_operator',s);
    verifyEqual(tc,C.dimension,n*q);
    verifyEqual(tc,C.steps,s);
    verifyTrue(tc,contains(C.approximation,'not_exact_expm'));
    F=C.apply;
    verifyEqual(tc,F([],'size'),[n*q n*q]);
    W=F(eye(n*q),'notransp');
    verifyLessThan(tc,norm(W'*W-eye(n*q),'fro'),4e-11);
    X=randn(n*q,3)+1i*randn(n*q,3);
    verifyLessThan(tc,norm(F(X,'notransp')-W*X,'fro'),3e-11);
    verifyLessThan(tc,norm(F(X,'ctransp')-W'*X,'fro'),3e-11);
    verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),3e-11);
    errors(j)=norm(W-U,'fro');
end
verifyGreaterThan(tc,errors(1)/errors(2),3.1);
verifyGreaterThan(tc,errors(2)/errors(3),3.1);
verifyGreaterThan(tc,errors(3),1e-8);
end

function test_large_sparse_approximation_and_errors(tc)
n=100;q=2;
C=graph_connection_propagator(n,q,[],[],[],1.1,'cayley_operator',8);
X=randn(n*q,2)+1i*randn(n*q,2);
F=C.apply;
verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),3e-10);
verifyError(tc,@()F(ones(n*q-1,1),'notransp'), ...
    'um_graph_connection_cayley:Dimension');
verifyError(tc,@()F(ones(n*q,1),'wrong'), ...
    'um_graph_connection_cayley:Flag');
verifyError(tc,@()graph_connection_propagator(4,1,[1 2;2 1]), ...
    'graph_connection_propagator:Edges');
verifyError(tc,@()graph_connection_propagator(4,1,[1 1]), ...
    'graph_connection_propagator:Edges');
verifyError(tc,@()graph_connection_propagator(4,1,[1 2],-1), ...
    'graph_connection_propagator:Weights');
verifyError(tc,@()graph_connection_propagator(4,2,[1 2],1,2*eye(2)), ...
    'graph_connection_propagator:Links');
verifyError(tc,@()graph_connection_propagator(4,1,[],[],[],1,'operator'), ...
    'graph_connection_propagator:Mode');
end
