function tests=test_opuc_christoffel
tests=functiontests(localfunctions);
end

function test_positive_reweighting_preserves_nodes_but_not_cyclic_weights(tc)
for n=[1 2 3 7 11]
    z=exp(2i*pi*((0:n-1)'+.23)/n);
    w=(1:n)'.^1.4;w=w/sum(w);
    zeta=.4+.3i;
    [U,info]=opuc_christoffel(n,z,w,zeta);
    verifyEqual(tc,info.dimension,n);
    verifyEmpty(tc,info.removed_indices);
    verifyLessThan(tc,norm(info.weights- ...
        (w.*abs(z-zeta).^2)/sum(w.*abs(z-zeta).^2)),3e-14);
    verifyLessThan(tc,norm(full(U'*U-eye(n)),'fro'),4e-10*n);
    [a,b]=eig(full(U));lam=diag(b);
    for j=1:n
        [dist,k]=min(abs(lam-z(j)));
        verifyLessThan(tc,dist,5e-9);
        verifyLessThan(tc,abs(abs(a(1,k))^2-info.weights(j)),5e-8);
    end
    S=opuc_christoffel(n,z,w,zeta,'cmv','spectral');
    verifyEqual(tc,S.dimension,n);
    verifyEqual(tc,S.finite.kind,'exact');
    verifyLessThan(tc,abs(S.finite.moments(2)-sum(z.^2)/n),3e-12);
    verifyLessThan(tc,abs(S.cyclic.moments(2)- ...
        sum(info.weights.*z.^2)),3e-12);
end
end

function test_exact_support_deletion_and_order_reduction(tc)
n=8;
z=exp(2i*pi*((0:n-1)'+.17)/n);w=(1:n)';w=w/sum(w);
[U,info]=opuc_christoffel(n,z,w,z(3),'hessenberg');
verifyEqual(tc,size(U),[n-1 n-1]);
verifyEqual(tc,info.dimension,n-1);
verifyEqual(tc,info.removed_indices,3);
verifyLessThan(tc,norm(full(U'*U-eye(n-1)),'fro'),5e-10);
verifyLessThan(tc,abs(sum(info.weights)-1),5e-14);
for j=[1:2 4:n]
    verifyLessThan(tc,min(abs(eig(full(U))-z(j))),6e-9);
end
S=opuc_christoffel(n,z,w,z(3),'hessenberg','spectral');
verifyEqual(tc,S.dimension,n-1);
verifyLessThan(tc,abs(S.finite.moments(1)- ...
    sum(info.nodes)/(n-1)),3e-13);
verifyEqual(tc,S.christoffel.removed_indices,3);
end

function test_trivial_factor_and_conjugation(tc)
n=6;
z=exp(2i*pi*((0:n-1)'+.1)/n);w=(1:n)';w=w/sum(w);
[U,info]=opuc_christoffel(n,z,w,0);
V=inverse_opuc(n,z,w);
verifyLessThan(tc,norm(full(U-V),'fro'),3e-12);
verifyLessThan(tc,norm(info.weights-w),3e-14);
[Uc,ic]=opuc_christoffel(n,conj(z),w,.25-.5i);
[Ut,it]=opuc_christoffel(n,z,w,.25+.5i);
verifyLessThan(tc,norm(full(Uc-conj(Ut)),'fro'),4e-11);
verifyLessThan(tc,norm(ic.weights-it.weights),4e-13);
F=opuc_christoffel(n,z,w,.2+.15i,'cmv','operator');
M=opuc_christoffel(n,z,w,.2+.15i);
X=randn(n,3)+1i*randn(n,3);
verifyEqual(tc,F([],'size'),[n n]);
verifyLessThan(tc,norm(F(X,'notransp')-M*X,'fro'),4e-11);
verifyLessThan(tc,norm(F(X,'ctransp')-M'*X,'fro'),4e-11);
verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),4e-11);
end

function test_annihilation_and_bad_inputs(tc)
verifyError(tc,@()opuc_christoffel(1,1,1,1), ...
    'opuc_christoffel:ZeroMeasure');
verifyError(tc,@()opuc_christoffel(2,[1 -1],[1 1],Inf), ...
    'opuc_christoffel:Zeta');
verifyError(tc,@()opuc_christoffel(3,[1 1 -1],[1 1 1],.2), ...
    'inverse_opuc:Distinct');
verifyError(tc,@()opuc_christoffel(2,[1 -1],[0 1],.2), ...
    'inverse_opuc:Weights');
end
