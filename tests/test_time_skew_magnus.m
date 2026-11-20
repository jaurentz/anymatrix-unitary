function tests=test_time_skew_magnus
tests=functiontests(localfunctions);
end

function test_stage_orthogonality_and_adjoint(tc)
for n=[3 4 7 11]
    for method={'exp_midpoint','cayley_midpoint'}
        m=method{1};T=1.37;
        S=time_skew_magnus(n,T,6,2.3,.85,m,'data');
        verifyTrue(tc,issparse(S.A0));
        verifyTrue(tc,issparse(S.A1));
        verifyLessThan(tc,norm(S.A0+S.A0','fro'),2e-15);
        verifyLessThan(tc,norm(S.A1+S.A1','fro'),2e-15);
        verifyGreaterThan(tc,norm(S.A0*S.A1-S.A1*S.A0,'fro'),.01);
        U=time_skew_magnus(n,T,6,2.3,.85,m);
        verifyEqual(tc,size(U),[n n]);
        verifyTrue(tc,isreal(U));
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),7e-12*n);
        verifyLessThan(tc,abs(det(U)-1),6e-12*n);
        F=time_skew_magnus(n,T,6,2.3,.85,m,'operator');
        X=randn(n,3)+1i*randn(n,3);
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),7e-12*n);
        verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),7e-12*n);
        verifyLessThan(tc,norm(F(X,'transp')-U'*X,'fro'),7e-12*n);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
            7e-12*n);
        reference=eye(n);
        for j=1:6
            A=S.at(S.midpoints(j));
            if strcmp(m,'exp_midpoint')
                Q=expm(S.step_size*full(A));
            else
                C=.5*S.step_size*full(A);
                Q=(eye(n)-C)\(eye(n)+C);
            end
            reference=Q*reference;
        end
        verifyLessThan(tc,norm(reference-U,'fro'),7e-12*n);
        R=time_skew_magnus(n,T,6,2.3,.85,m,'spectral');
        verifyEqual(tc,R.finite.kind,'unknown');
        verifyEmpty(tc,R.finite.eigenvalues);
    end
end
end

function test_second_order_noncommuting_convergence(tc)
n=4;T=1.37;nu=2.3;amp=.85;
P=time_skew_magnus(n,T,4,nu,amp,'exp_midpoint','data');
reference=ode45(@(t,y)reshape(P.at(t)*reshape(y,n,n),n*n,1), ...
    [0 T],reshape(eye(n),n*n,1), ...
    odeset('RelTol',1e-12,'AbsTol',1e-13));
Y=reshape(deval(reference,T),n,n);
for method={'exp_midpoint','cayley_midpoint'}
    m=method{1};errors=zeros(1,3);
    for k=1:3
        steps=4*2^(k-1);
        U=time_skew_magnus(n,T,steps,nu,amp,m);
        errors(k)=norm(U-Y,'fro');
    end
    verifyGreaterThan(tc,errors(1)/errors(2),3.5);
    verifyGreaterThan(tc,errors(2)/errors(3),3.5);
    verifyLessThan(tc,errors(3),.002);
end
end

function test_constant_spectral_and_not_generically_symplectic(tc)
n=4;J=[zeros(2),eye(2);-eye(2),zeros(2)];
U=time_skew_magnus(n,1.37,8,2.3,.85,'cayley_midpoint');
verifyGreaterThan(tc,norm(U'*J*U-J,'fro'),.05);
for method={'exp_midpoint','cayley_midpoint'}
    m=method{1};
    n=7;T=1.21;
    U=time_skew_magnus(n,T,9,2.3,0,m);
    S=time_skew_magnus(n,T,9,2.3,0,m,'spectral');
    verifyEqual(tc,S.finite.kind,'exact');
    verifyEqual(tc,S.dimension,n);
    z=S.finite.eigenvalues();
    for a=eig(U).'
        verifyLessThan(tc,min(abs(a-z)),3e-11);
    end
    for k=-5:5
        verifyLessThan(tc,abs(S.finite.moments(k)-trace(U^k)/n), ...
            3e-11);
    end
    verifyError(tc,@()S.finite.moments(.5), ...
        'um_skew_moments:Moments');
end
end

function test_large_sparse_operator_and_input(tc)
n=128;
F=time_skew_magnus(n,1.1,5,2.1,.8,'cayley_midpoint','operator');
X=randn(n,2)+1i*randn(n,2);
verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
    1e-10);
verifyError(tc,@()F(ones(n-1,2),'notransp'), ...
    'um_skew_midpoint_apply:Dimension');
verifyError(tc,@()F(ones(n,2),'invalid'), ...
    'um_skew_midpoint_apply:Flag');
verifyError(tc,@()time_skew_magnus(4,1,8,1,1,'bad'), ...
    'time_skew_magnus:Method');
verifyError(tc,@()time_skew_magnus(4,1,8,1,1,[], 'bad'), ...
    'time_skew_magnus:Mode');
verifyError(tc,@()time_skew_magnus(129,1,8,1,1,'exp_midpoint'), ...
    'time_skew_magnus:Order');
end
