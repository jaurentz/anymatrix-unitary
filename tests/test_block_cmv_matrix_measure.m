function tests=test_block_cmv_matrix_measure
tests=functiontests(localfunctions);
end

function test_matrix_free_noncommuting_cores_and_pencil(tc)
for q=[1 2 3]
    for n=[1 2 4 7]
        alpha=complex(zeros(q,q,max(n-1,0)));
        for k=1:n-1
            D=diag(exp(1i*(1:q)*(.13+k*.17)));
            J=fliplr(eye(q));
            alpha(:,:,k)=.11*D+.09i*J*D;
        end
        beta=diag(exp(1i*(1:q)*.27));
        U=block_cmv_pencil(n,q,alpha,beta);
        P=block_cmv_pencil(n,q,alpha,beta,'pencil');
        F=block_cmv_pencil(n,q,alpha,beta,'operator');
        X=randn(n*q,4)+1i*randn(n*q,4);
        verifyTrue(tc,issparse(U));
        verifyEqual(tc,F([],'size'),[n*q n*q]);
        verifyLessThan(tc,norm(full(U'*U-eye(n*q)),'fro'),2e-10*n*q);
        verifyLessThan(tc,norm(full(P.B\P.A-U),'fro'),2e-10*n*q);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),2e-10*n*q);
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),2e-10*n*q);
        verifyLessThan(tc,norm(F(X,'transp')-U'*X,'fro'),2e-10*n*q);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
            2e-10*n*q);
        if q==1
            sc=[reshape(alpha,[],1);beta];
            C=finite_cmv(n,sc);
            verifyLessThan(tc,norm(full(U-C),'fro'),2e-10*n);
        end
    end
end
end

function test_matrix_valued_quadrature_vs_scalar_trace(tc)
n=5;q=2;
alpha=zeros(q,q,n-1);
alpha(:,:,1)=[.18, .12i;-.05,.09];
alpha(:,:,2)=[-.04,.21; .1i,.03];
alpha(:,:,3)=[.02i,.15; -.11,.07];
alpha(:,:,4)=[.08,-.03; .18i,-.02];
verifyGreaterThan(tc,norm(alpha(:,:,1)*alpha(:,:,2)- ...
    alpha(:,:,2)*alpha(:,:,1),'fro'),.005);
beta=[0 1;-1 0];
U=full(block_cmv_pencil(n,q,alpha,beta));
Q=block_cmv_pencil(n,q,alpha,beta,'quadrature');
verifyEqual(tc,Q.kind,'numerical');
verifyEqual(tc,Q.dimension,n*q);
verifyEqual(tc,size(Q.weights),[q q n*q]);
verifyLessThan(tc,norm(sum(Q.weights,3)-eye(q),'fro'),1e-10);
verifyLessThan(tc,abs(sum(Q.trace_weights)-1),2e-14);
E=zeros(n*q,q);E(1:q,:)=eye(q);
for j=1:n*q
    W=Q.weights(:,:,j);
    verifyLessThan(tc,norm(W-W','fro'),3e-12);
    verifyGreaterThanOrEqual(tc,min(eig(W)),-3e-12);
end
ks=[-4 -2 -1 0 1 2 5];
Mc=Q.cyclic_moments(ks);
Mt=Q.trace_moments(ks);
for j=1:numel(ks)
    k=ks(j);
    verifyLessThan(tc,norm(Mc(:,:,j)-E'*U^k*E,'fro'), ...
        4e-10);
    verifyLessThan(tc,abs(Mt(j)-trace(U^k)/(n*q)),4e-10);
end
verifyGreaterThan(tc,norm(Mc(:,:,2)-Mt(2)*eye(q),'fro'),1e-3);
S=block_cmv_pencil(n,q,alpha,beta,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyError(tc,@()Q.cyclic_moments(.5), ...
    'um_block_cmv_moments:Order');
end

function test_large_matrix_free_and_invalid(tc)
n=120;q=2;
F=block_cmv_pencil(n,q,[],[],'operator');
X=randn(n*q,3)+1i*randn(n*q,3);
verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
    1e-10);
verifyError(tc,@()F(ones(n*q-1,1),'notransp'), ...
    'um_block_cmv_apply:Dimension');
verifyError(tc,@()F(ones(n*q,1),'bad'), ...
    'um_block_cmv_apply:Flag');
verifyError(tc,@()block_cmv_pencil(70,2,[],[],'quadrature'), ...
    'block_cmv_pencil:QuadratureSize');
verifyError(tc,@()block_cmv_pencil(3,2,[],[],'bad'), ...
    'block_cmv_pencil:Mode');
verifyError(tc,@()block_cmv_pencil(3,2,repmat(2*eye(2),1,1,2),eye(2)), ...
    'block_cmv_pencil:Alpha');
end
