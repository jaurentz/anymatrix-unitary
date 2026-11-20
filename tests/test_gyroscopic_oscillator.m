function tests=test_gyroscopic_oscillator
%TEST_GYROSCOPIC_OSCILLATOR Gyro-QEP, metric and canonical identities.
tests=functiontests(localfunctions);
end

function test_coefficients_and_invariants(tc)
for n=[2 4 7]
    for method={'midpoint','gauss2','gauss3'}
        s=method{1};
        P=gyroscopic_oscillator(n,.3,1.1,.7,.4,s,'pencil');
        q=gyroscopic_oscillator(n,.3,1.1,.7,.4,s,'quadratic');
        M=q.A2;K=q.A0;C=q.A1;
        verifyLessThan(tc,norm(M-M','fro'),1e-13);
        verifyLessThan(tc,norm(K-K','fro'),1e-13);
        verifyLessThan(tc,norm(C+C','fro'),1e-13);
        verifyGreaterThan(tc,norm(C,'fro'),.01);
        verifyGreaterThan(tc,min(eig(full(M))),0);
        verifyGreaterThan(tc,min(eig(full(K))),0);
        verifyGreaterThan(tc,norm(M*K-K*M,'fro'),1e-4);
        G=P.G;S=P.Mstep;T=P.T;
        verifyLessThan(tc,norm(P.physicalGenerator'*G+G*P.physicalGenerator,'fro'),1e-11);
        verifyLessThan(tc,norm(S'*G*S-G,'fro'),5e-10*norm(G,'fro'));
        verifyLessThan(tc,norm(T'*T-G,'fro'),2e-11*norm(G,'fro'));
        verifyLessThan(tc,norm(P.U'*P.U-eye(2*n),'fro'),5e-10*n);
        verifyLessThan(tc,norm(T*S/T-P.U,'fro'),5e-10*n);
        verifyLessThan(tc,norm(P.A-G*S,'fro'),2e-11*norm(G,'fro'));
        verifyLessThan(tc,norm(P.B-G,'fro'),2e-11*norm(G,'fro'));
        verifyLessThan(tc,norm(P.A'*(P.B\P.A)-P.B,'fro'),1e-9*norm(G,'fro'));
        J=[zeros(n),eye(n);-eye(n),zeros(n)];
        B=P.canonicalMap;Y=B*S/B;
        verifyLessThan(tc,norm(Y'*J*Y-J,'fro'),5e-10*n);
        Jv=B'*J*B;
        verifyLessThan(tc,norm(S'*Jv*S-Jv,'fro'),5e-10*n);
        verifyEqual(tc,P.classification,'circle_spectrum_only');
        lam=eig(full(P.A),full(P.B));
        verifyLessThan(tc,max(abs(abs(lam)-1)),1e-8);
    end
end
end

function test_continuous_qep_and_discrete_pade(tc)
n=6;h=.35;
q=gyroscopic_oscillator(n,h,1.3,.8,.5,'gauss2','quadratic');
P=gyroscopic_oscillator(n,h,1.3,.8,.5,'gauss2','pencil');
lam=polyeig(full(q.A0),full(q.A1),full(q.A2));
eA=eig(full(P.physicalGenerator));
eU=eig(P.U);
verifyEqual(tc,numel(lam),2*n);
verifyLessThan(tc,max(abs(real(lam))),3e-9);
for k=1:numel(lam)
    z=lam(k);
    verifyLessThan(tc,min(abs(eA-z)),1e-8);
    rz=(1+h*z/2+(h*z)^2/12)/(1-h*z/2+(h*z)^2/12);
    verifyLessThan(tc,min(abs(eU-rz)),1e-8);
end
end

function test_actions_and_zero_gyro_limit(tc)
for gamma=[0 .4]
    for method={'midpoint','gauss2','gauss3'}
        U=gyroscopic_oscillator(8,.3,1,.5,gamma,method{1});
        O=gyroscopic_oscillator(8,.3,1,.5,gamma,method{1},'operator');
        X=randn(16,3)+1i*randn(16,3);
        verifyEqual(tc,O([],'size'),[16 16]);
        verifyLessThan(tc,norm(O(X,'notransp')-U*X,'fro'),1e-9);
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'),1e-9);
        verifyLessThan(tc,norm(O(X,'transp')-U.'*X,'fro'),1e-9);
        verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X,'fro'),1e-9);
    end
end
end

function test_properties_unknown_spectrum_and_rejection(tc)
[A,props]=gyroscopic_oscillator();
verifyEmpty(tc,A);
verifyTrue(tc,any(strcmp(props,'unitary')));
verifyTrue(tc,any(strcmp(props,'orthogonal')));
v=gyroscopic_oscillator(4,[],[],[],[],[],'spectral');
verifyEqual(tc,v.finite.kind,'unknown');
verifyEqual(tc,v.dimension,8);
verifyError(tc,@()gyroscopic_oscillator(4,.25,1,1,.4,'bad'), ...
    'gyroscopic_oscillator:Method');
verifyError(tc,@()gyroscopic_oscillator(4,.25,1,1,.4,'midpoint','bad'), ...
    'gyroscopic_oscillator:Mode');
verifyError(tc,@()gyroscopic_oscillator(4,.25,-1), ...
    'MATLAB:expectedPositive');
O=gyroscopic_oscillator(4,[],[],[],[],[],'operator');
verifyError(tc,@()O(ones(7,1),'notransp'), ...
    'um_gyroscopic_oscillator_apply:Dimension');
verifyError(tc,@()O(ones(8,1),'invalid'), ...
    'um_gyroscopic_oscillator_apply:Flag');
end
