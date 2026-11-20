function tests=test_hahn_dual_geodesic
%TEST_HAHN_DUAL_GEODESIC Independent Hahn-duality and geodesic identities.
tests=functiontests(localfunctions);
end

function test_hahn_orthogonality_jacobi_and_duality(testCase)
for n=[1 2 3 4 7 12 20]
    for ab=[0 .4 -0.75;.3 .9 .2]
        a=ab(1);b=ab(2);
        H=hahn_transform(n,a,b);
        D=dual_hahn_transform(n,a,b);
        verifyEqual(testCase,size(H),[n n]);
        verifyTrue(testCase,isreal(H));
        verifyLessThan(testCase,norm(H'*H-eye(n),'fro'),3e-11*n);
        verifyLessThan(testCase,norm(D'*D-eye(n),'fro'),3e-11*n);
        verifyLessThan(testCase,norm(D-H','fro'),1e-14*n);
        verifyGreaterThan(testCase,min(H(:,1)),0);
        [V,J,lambda]=um_hahn_basis(n,a,b);
        expected=(0:n-1)'.*((0:n-1)'+a+b+1);
        verifyLessThan(testCase,norm(lambda-expected),2e-11*n*n);
        verifyLessThan(testCase,norm(J*V-V*diag(expected),'fro'),2e-11*n*n);
        verifyLessThan(testCase,norm(H-V','fro'),1e-14*n);
        S=hahn_transform(n,a,b,'spectral');
        T=dual_hahn_transform(n,a,b,'spectral');
        verifyEqual(testCase,S.finite.kind,'unknown');
        verifyEqual(testCase,T.finite.kind,'unknown');
    end
end
end

function test_hahn_explicit_polynomial_and_weights(testCase)
for n=[3 5 8 13]
    N=n-1; x=(0:N)';
    for a=[0 .4 1.2]
        for b=[0 .8 1.5]
            H=hahn_transform(n,a,b);
            logw=gammaln(a+x+1)-gammaln(a+1)-gammaln(x+1) ...
                +gammaln(b+N-x+1)-gammaln(b+1)-gammaln(N-x+1);
            w=exp(logw-max(logw));
            v0=sqrt(w/sum(w));
            verifyLessThan(testCase,norm(H(1,:)'-v0),6e-11*n);
            poly1=1-(a+b+2)*x/((a+1)*N);
            observed=(H(2,:) ./ H(1,:))';
            observed=observed/observed(1);
            verifyLessThan(testCase,norm(observed-poly1),2e-9*n);
        end
    end
end
end

function test_hahn_random_parameter_eigenphase_conjugacy(testCase)
for n=[2 3 6 9 15]
    H=hahn_transform(n,.35,1.2);
    D=dual_hahn_transform(n,.35,1.2);
    z=eig(H);
    w=eig(D);
    verifyLessThan(testCase,max(arrayfun(@(a)min(abs(a-conj(w))),z)),2e-10*n);
    verifyLessThan(testCase,max(abs(abs(z)-1)),5e-11*n);
    verifyLessThan(testCase,max(abs(abs(w)-1)),5e-11*n);
end
end

function test_hahn_validation_and_order_one(testCase)
verifyEqual(testCase,hahn_transform(1),1);
verifyEqual(testCase,dual_hahn_transform(1),1);
verifyError(testCase,@()hahn_transform(5,-1,.2),'MATLAB:notGreater');
verifyError(testCase,@()dual_hahn_transform(5,.2,-1),'MATLAB:notGreater');
verifyError(testCase,@()hahn_transform(5,.2,.4,'operator'),'hahn_transform:Mode');
verifyError(testCase,@()dual_hahn_transform(5,.2,.4,'operator'),'dual_hahn_transform:Mode');
end

function test_geodesic_rotation_orthogonality_endpoints(testCase)
for n=[2 3 4 8 17]
    for theta=[0 .17 .7 1.3 2.5]
        x=zeros(n,1);x(1)=1;
        y=zeros(n,1);y(1)=cos(theta);y(2)=sin(theta);
        U=sphere_geodesic_transport(n,x,y);
        verifyTrue(testCase,isreal(U));
        verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),3e-11*n);
        verifyLessThan(testCase,norm(U*x-y),3e-11*n);
        verifyLessThan(testCase,abs(det(U)-1),3e-10*n);
        verifyLessThan(testCase,norm(U(3:end,3:end)-eye(n-2),'fro'),3e-11*n);
        % The restriction to tangent vectors is true Levi-Civita transport.
        v=zeros(n,1);v(2)=1;
        expected=v-(y'*v)/(1+x'*y)*(x+y);
        verifyLessThan(testCase,norm(U*v-expected),4e-11*n);
        verifyLessThan(testCase,abs(y'*(U*v)),3e-11*n);
        S=sphere_geodesic_transport(n,x,y,'spectral');
        verifyEqual(testCase,S.finite.kind,'exact');
        z=S.finite.eigenvalues();
        actual=eig(U);
        verifyLessThan(testCase,max(arrayfun(@(q)min(abs(q-z)),actual)),3e-10*n);
        for k=-3:5
            verifyLessThan(testCase,abs(S.finite.moments(k)-trace(U^k)/n),4e-10*n);
        end
        verifyLessThan(testCase,abs(sum(S.finite.atoms.weights)-1),1e-14);
    end
end
end

function test_geodesic_general_endpoints_and_operator(testCase)
for n=[2 3 5 12 30]
    x=(1:n)';x=x/norm(x);
    z=sin((1:n)'*.83);
    z=z-x*(x'*z);z=z/norm(z);
    y=cos(.65)*x+sin(.65)*z;
    U=sphere_geodesic_transport(n,x,y);
    F=sphere_geodesic_transport(n,x,y,'operator');
    verifyEqual(testCase,F([],'size'),[n n]);
    X=randn(n,3)+1i*randn(n,3);
    verifyLessThan(testCase,norm(U'*U-eye(n),'fro'),2e-10*n);
    verifyLessThan(testCase,norm(U*x-y),2e-10*n);
    verifyLessThan(testCase,norm(F(X,'notransp')-U*X,'fro'),2e-10*n);
    verifyLessThan(testCase,norm(F(X,'ctransp')-U'*X,'fro'),2e-10*n);
    verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),2e-10*n);
end
end

function test_geodesic_large_implicit_and_degenerate(testCase)
n=10000;
F=sphere_geodesic_transport(n,[],[],'operator');
X=randn(n,2)+1i*randn(n,2);
verifyLessThan(testCase,norm(F(F(X,'notransp'),'transp')-X,'fro'),3e-9);
S=sphere_geodesic_transport(n,[],[],'spectral');
verifyEqual(testCase,S.dimension,n);
verifyEqual(testCase,numel(S.finite.eigenvalues()),n);
x=zeros(4,1);x(1)=1;
I=sphere_geodesic_transport(4,x,x);
verifyLessThan(testCase,norm(I-eye(4),'fro'),1e-14);
SI=sphere_geodesic_transport(4,x,x,'spectral');
verifyEqual(testCase,SI.finite.atoms.weights,1);
end

function test_geodesic_errors(testCase)
verifyError(testCase,@()sphere_geodesic_transport(1),'MATLAB:notGreaterEqual');
verifyError(testCase,@()sphere_geodesic_transport(4,[1;0;0],[]), ...
    'sphere_geodesic_transport:Start');
verifyError(testCase,@()sphere_geodesic_transport(4,[2;0;0;0],[]), ...
    'sphere_geodesic_transport:Start');
x=[1;0;0];y=-x;
verifyError(testCase,@()sphere_geodesic_transport(3,x,y), ...
    'sphere_geodesic_transport:Antipodal');
verifyError(testCase,@()sphere_geodesic_transport(3,x,[], 'bad'), ...
    'sphere_geodesic_transport:Mode');
F=sphere_geodesic_transport(3,[],[],'operator');
verifyError(testCase,@()F(ones(2,1),'notransp'),'um_geodesic_apply:Dimension');
verifyError(testCase,@()F(ones(3,1),'bad'),'um_geodesic_apply:Flag');
end
