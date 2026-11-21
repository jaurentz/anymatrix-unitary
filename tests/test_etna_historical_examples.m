function tests=test_etna_historical_examples
%TEST_ETNA_HISTORICAL_EXAMPLES ETNA 44 (2015), Section 7.1 and 7.5.
tests=functiontests(localfunctions);
end

function test_example71_correlated_core_columns(tc)
for n=[2 3 8 19]
    for B=[1 2 2.5]
        P=etna_hessenberg(n,'7.1',19,B,'parameters');
        m=P.common_m;p=P.common_p;
        k=(1:n-1)';
        col1=(-1).^(k-1).*p.*sqrt(1-m.^(2*B.^(-k)));
        col2=m.^(B.^(-k));
        verifyGreaterThan(tc,m,0);
        verifyLessThan(tc,m,1);
        verifyLessThan(tc,abs(abs(p)-1),3e-15);
        verifyLessThan(tc,norm(P.first_core_columns(:,1)-col1), ...
            3e-14*n);
        verifyLessThan(tc,norm(P.first_core_columns(:,2)-col2), ...
            3e-14*n);
        verifyLessThan(tc,norm(P.alpha-conj(col1)),3e-14*n);
        verifyEqual(tc,P.beta,(-1)^(n-1));
        U=full(etna_hessenberg(n,'7.1',19,B));
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),4e-11*n);
        verifyLessThan(tc,norm(tril(U,-2),'fro'),2e-11*n);
        F=etna_hessenberg(n,'7.1',19,B,'operator');
        X=randn(n,3)+1i*randn(n,3);
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),3e-11*n);
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),3e-11*n);
    end
end
end

function test_example71_real_and_seeded_replay(tc)
for n=[3 8 20]
    A=full(etna_hessenberg(n,'7.1-real',5,2));
    B=full(etna_hessenberg(n,'7.1-real',5,2));
    verifyEqual(tc,A,B,'AbsTol',0);
    verifyTrue(tc,isreal(A));
    verifyLessThan(tc,norm(A'*A-eye(n),'fro'),3e-11*n);
    P=etna_hessenberg(n,'7.1-real',5,2,'parameters');
    verifyEqual(tc,P.common_p,1);
    verifyLessThan(tc,max(abs(imag(P.alpha))),1e-14);
end
old=rng;
restore=onCleanup(@()rng(old)); %#ok<NASGU>
rng(38);before=rng;
etna_hessenberg(13,'7.1',33,2);
after=rng;
verifyEqual(tc,before,after);
end

function test_inverse_example75_prescribed_circle_spectrum(tc)
for n=[4 8 12 20 28]
    for kind={'7.5','7.5-real'}
        name=kind{1};seed=11;
        S=etna_hessenberg(n,name,seed,[],'spectral');
        P=etna_hessenberg(n,name,seed,[],'parameters');
        U=full(etna_hessenberg(n,name,seed));
        nodes=exp(2i*pi*((0:n-1)'+.5)/n);
        verifyEqual(tc,S.family,'etna_hessenberg');
        verifyEqual(tc,S.finite.kind,'exact');
        verifyLessThan(tc,norm(S.finite.eigenvalues()-nodes),3e-13*n);
        verifyLessThan(tc,abs(sum(P.weights)-1),1e-12);
        verifyGreaterThan(tc,min(P.weights),0);
        verifyLessThan(tc,abs(P.beta+1),1e-11*n);
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),5e-10*n);
        verifyLessThan(tc,norm(tril(U,-2),'fro'),5e-10*n);
        lam=eig(U);
        for j=1:n
            verifyLessThan(tc,min(abs(lam-nodes(j))),3e-8*n);
        end
        verifyLessThan(tc,P.maximum_forward_node_error,3e-8*n);
        verifyTrue(tc,P.eigenvalue_matching_unique);
        verifyLessThan(tc,P.maximum_forward_weight_error,2e-7*n);
        % Independent diagonalization: compare cyclic eigenvector weights,
        % NOT the equal-weight empirical trace eigenphase measure.
        [V,D]=eig(U);
        eigenvalues=diag(D);
        distance=abs(nodes-eigenvalues.');
        [~,match]=min(distance,[],2);
        verifyEqual(tc,numel(unique(match)),n);
        numericalCyclicWeights=abs(V(1,match).').^2;
        verifyLessThan(tc,max(abs(numericalCyclicWeights-P.weights)), ...
            2e-7*n);
        for k=0:3
            observedMoment=sum(numericalCyclicWeights.* ...
                eigenvalues(match).^k);
            targetMoment=sum(P.weights.*nodes.^k);
            verifyLessThan(tc,abs(observedMoment-targetMoment), ...
                3e-7*n);
        end
        F=etna_hessenberg(n,name,seed,[],'operator');
        X=randn(n,2)+1i*randn(n,2);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'), ...
            2e-10*n);
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'), ...
            2e-10*n);
        if strcmp(name,'7.5-real')
            verifyLessThan(tc,norm(imag(U),'fro'),1e-8*n);
        end
    end
end
end

function test_historical_and_modern_example75_are_distinct(tc)
A=etna_hessenberg(10,'7.5',2,[],'parameters');
B=etna_hessenberg(10,'7.5',5,[],'parameters');
verifyGreaterThan(tc,norm(A.weights-B.weights),.01);
verifyEqual(tc,A.paper_example,'7.5-modern-analogue');
verifyTrue(tc,A.eigenvalue_matching_unique);
verifyLessThan(tc,A.maximum_forward_weight_error,2e-7*10);
verifyEqual(tc,etna_hessenberg(10,'7.5',2,[],'spectral').finite.kind,'exact');
S=etna_hessenberg(10,'7.1',2,2,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyNotEmpty(tc,S.notes);
verifyError(tc,@()etna_hessenberg(5,'7.1',2,.5), ...
    'etna_hessenberg:BetaPower');
verifyError(tc,@()etna_hessenberg(5,'new',2), ...
    'etna_hessenberg:Kind');
end
