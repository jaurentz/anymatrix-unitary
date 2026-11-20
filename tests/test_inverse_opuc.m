function tests=test_inverse_opuc
tests=functiontests(localfunctions);
end

function test_nodes_weights_unitarity_and_multiple_representations(tc)
for n=[1 2 3 5 8 12]
    theta=(0:n-1)'*2*pi/n+.1+.07*sin((1:n)');
    z=exp(1i*theta);
    weights=(1:n)'.^1.5;
    weights=weights/sum(weights);
    params=inverse_opuc(n,z,weights,'cmv','parameters');
    verifyEqual(tc,numel(params.alpha),n-1);
    verifyLessThan(tc,abs(params.beta- ...
        (-1)^(n+1)*conj(prod(z))),2e-12);
    verifyLessThan(tc,max([0;abs(params.alpha)]) ,1);
    for rep={'cmv','hessenberg','inverse_hessenberg'}
        name=rep{1};
        U=full(inverse_opuc(n,z,weights,name));
        verifyEqual(tc,size(U),[n n]);
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),5e-11*n);
        [a,b]=eig(U);lam=diag(b);
        for k=1:n
            [dist,j]=min(abs(lam-z(k)));
            verifyLessThan(tc,dist,3e-9);
            verifyLessThan(tc,abs(abs(a(1,j))^2-weights(k)),3e-8);
        end
        if strcmp(name,'cmv')
            verifyLessThanOrEqual(tc,nnz(U),5*n);
        elseif strcmp(name,'hessenberg')
            verifyLessThan(tc,norm(tril(U,-2),'fro'),3e-12*n);
        end
        F=inverse_opuc(n,z,weights,name,'operator');
        X=randn(n,3)+1i*randn(n,3);
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),4e-11*n);
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),4e-11*n);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
            5e-11*n);
    end
end
end

function test_trace_vs_cyclic_moments(tc)
n=7;
z=exp(2i*pi*((0:n-1)'+.23)/n);
w=(1:n)';w=w/sum(w);
S=inverse_opuc(n,z,w,'cmv','spectral');
verifyEqual(tc,S.finite.kind,'exact');
verifyEqual(tc,S.dimension,n);
verifyEqual(tc,S.cyclic.kind,'exact');
verifyEqual(tc,numel(S.finite.eigenvalues()),n);
verifyLessThan(tc,norm(S.cyclic.weights-w),1e-13);
for k=-5:5
    verifyLessThan(tc,abs(S.finite.moments(k)-sum(z.^k)/n),3e-14);
    verifyLessThan(tc,abs(S.cyclic.moments(k)-sum(w.*z.^k)),3e-14);
end
verifyGreaterThan(tc,abs(S.cyclic.moments(1)-S.finite.moments(1)),.03);
params=inverse_opuc(n,z,w,'cmv','parameters');
verifyLessThan(tc,abs(params.alpha(1)-conj(sum(w.*z))),1e-12);
verifyError(tc,@()S.finite.moments(.5),'inverse_opuc:Moments');
end

function test_clustered_and_conjugation_cases(tc)
n=6;
theta=[-.1; -.1+.008; .7;1.9;2.7;4.2];
z=exp(1i*theta);
w=[1;2;3;3;2;1]/12;
U=inverse_opuc(n,z,w);
verifyLessThan(tc,norm(full(U'*U-eye(n)),'fro'),1e-9);
for j=1:n
    verifyLessThan(tc,min(abs(eig(full(U))-z(j))),1e-8);
end
Uc=inverse_opuc(n,conj(z),w);
verifyLessThan(tc,norm(full(Uc-conj(U)),'fro'),5e-9);
end

function test_defaults_and_bad_spectral_data(tc)
[U,p]=inverse_opuc();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyEqual(tc,size(inverse_opuc(5)),[5 5]);
verifyError(tc,@()inverse_opuc(3,[1 1 -1]),'inverse_opuc:Distinct');
verifyError(tc,@()inverse_opuc(3,[1 2 -1]),'inverse_opuc:Nodes');
verifyError(tc,@()inverse_opuc(2,[1 -1],[0 1]), ...
    'inverse_opuc:Weights');
verifyError(tc,@()inverse_opuc(2,[1 -1],[1 -2]), ...
    'inverse_opuc:Weights');
verifyError(tc,@()inverse_opuc(2,[1 -1],[1 1],'bad'), ...
    'inverse_opuc:Representation');
verifyError(tc,@()inverse_opuc(2,[1 -1],[1 1],'cmv','bad'), ...
    'inverse_opuc:Mode');
end
