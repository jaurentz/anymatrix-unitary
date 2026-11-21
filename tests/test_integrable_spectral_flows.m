function tests=test_integrable_spectral_flows
%TEST_INTEGRABLE_SPECTRAL_FLOWS Independent finite spectral and QR checks.
tests=functiontests(localfunctions);
end

function test_discrete_schur_isospectral_and_cyclic(tc)
n=7;z=exp(2i*pi*((0:n-1)'+.13)/n);
w=(1:n)';w=w/sum(w);h=.08;
for k=[0 1 3 16]
    S=discrete_schur_cmv(n,k,h,z,w,'spectral');
    g=1+2*h*real(z);
    target=w.*g.^k;target=target/sum(target);
    verifyLessThan(tc,norm(S.cyclic.weights-target),2e-13);
    verifyLessThan(tc,abs(S.finite.moments(-2)-sum(z.^(-2))/n),2e-13);
    verifyLessThan(tc,abs(S.cyclic.moments(2)-sum(target.*z.^2)),2e-13);
    verifyEqual(tc,S.finite.kind,'exact');
    U=full(discrete_schur_cmv(n,k,h,z,w));
    verifyLessThan(tc,norm(U'*U-eye(n),'fro'),2e-9);
    lam=eig(U);
    for j=1:n
        verifyLessThan(tc,min(abs(lam-z(j))),3e-8);
    end
    par=discrete_schur_cmv(n,k,h,z,w,'parameters');
    verifyLessThan(tc,norm(par.weights-target),2e-13);
    verifyLessThan(tc,abs(abs(par.beta)-1),1e-12);
end
end

function test_discrete_schur_limit_real_and_operator(tc)
z=[1;exp(.7i);exp(-.7i);-1];w=[.2;.3;.3;.2];n=4;
t=.5;err=zeros(3,1);
for j=1:3
    h=.2/2^(j-1);k=round(t/h);
    % Use t=0.4 to make all K*H exactly equal; compare at fixed t.
    h=.1/2^(j-1);k=round(.4/h);
    S=discrete_schur_cmv(n,k,h,z,w,'spectral');
    target=w.*exp(.8*real(z));target=target/sum(target);
    err(j)=norm(S.cyclic.weights-target);
end
verifyLessThan(tc,err(3),err(2));
verifyLessThan(tc,err(2),err(1));
A=full(discrete_schur_cmv(n,4,.1,z,w));
verifyLessThan(tc,norm(imag(A),'fro'),2e-10);
F=discrete_schur_cmv(n,4,.1,z,w,'operator');
X=randn(n,2)+1i*randn(n,2);
verifyEqual(tc,F([],'size'),[n n]);
verifyLessThan(tc,norm(F(X,'notransp')-A*X,'fro'),2e-10);
verifyLessThan(tc,norm(F(X,'ctransp')-A'*X,'fro'),2e-10);
S=discrete_schur_cmv(4,180,.1,[1;1i;-1;-1i],w,'spectral');
verifyGreaterThan(tc,S.cyclic.weights(1),.999);
verifyError(tc,@()discrete_schur_cmv(n,1,.5,z,w), ...
    'MATLAB:expectedLessThan'); % validateattributes bound
end

function test_toda_jacobi_qr_and_cayley(tc)
n=7;
for tau=[0 .3 -.4]
    P=toda_cayley(n,tau,.7,.2,.4,'parameters');
    J=full(P.jacobi);
    J0=full(P.initial_jacobi);
    verifyLessThan(tc,norm(J-J.','fro'),1e-12);
    verifyLessThan(tc,norm(tril(J,-2),'fro'),2e-12);
    verifyGreaterThan(tc,min(diag(J,-1)),0);
    verifyLessThan(tc,norm(J-P.orthogonal_factor'*J0* ...
        P.orthogonal_factor,'fro'),1e-10);
    verifyLessThan(tc,norm(sort(eig(J))-sort(eig(J0))),3e-10);
    U=toda_cayley(n,tau,.7,.2,.4);
    verifyLessThan(tc,norm(U'*U-eye(n),'fro'),2e-10);
    verifyLessThan(tc,norm(U-U.','fro'),2e-10);
    verifyLessThan(tc,norm(sort(angle(eig(U)))- ...
        sort(angle((eig(J0)-.7i)./(eig(J0)+.7i)))),3e-8);
    S=toda_cayley(n,tau,.7,.2,.4,'spectral');
    target=P.initial_weights.*exp(2*tau*P.jacobi_eigenvalues);
    target=target/sum(target);
    verifyLessThan(tc,norm(S.cyclic.weights-target),3e-12);
    verifyLessThan(tc,abs(S.finite.moments(1)-trace(U)/n),3e-11);
    for j=1:n
        [V,D]=eig(J);[~,id]=min(abs(diag(D)-P.jacobi_eigenvalues(j)));
        verifyLessThan(tc,abs(abs(V(1,id))^2-S.cyclic.weights(j)),2e-8);
    end
    F=toda_cayley(n,tau,.7,.2,.4,'operator');
    X=randn(n,2)+1i*randn(n,2);
    verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),3e-10);
    verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),3e-10);
end
end

function test_toda_custom_recurrence_and_validation(tc)
n=4;rec=struct('diagonal',[-1;-.2;.3;1], ...
    'subdiagonal',[.1;.4;.2]);
J=toda_cayley(n,.3,1,rec,[],'jacobi');
verifyLessThan(tc,norm(J-J','fro'),1e-13);
P=toda_cayley(n,.3,1,rec,[],'parameters');
verifyLessThan(tc,norm(J-P.jacobi,'fro'),1e-13);
verifyError(tc,@()toda_cayley(3,0,-.2),'MATLAB:expectedPositive');
verifyError(tc,@()toda_cayley(4,0,1, ...
    struct('diagonal',ones(4,1),'subdiagonal',[1;0;1])), ...
    'toda_cayley:Recurrence');
end

function test_hessenberg_qr_vs_inverse_cyclic_spectrum(tc)
n=6;z=exp(2i*pi*((0:n-1)'+.1)/n);
w=(1:n)';w=w/sum(w);
for t=[0 .3 -.5]
    H=schur_hessenberg_flow(n,t,z,w);
    I=schur_hessenberg_flow(n,t,z,w,'hessenberg','inverse');
    verifyLessThan(tc,norm(H'*H-eye(n),'fro'),2e-10);
    verifyLessThan(tc,norm(tril(H,-2),'fro'),1e-10);
    verifyLessThan(tc,norm(H-I,'fro'),2e-8);
    S=schur_hessenberg_flow(n,t,z,w,'hessenberg','qr','spectral');
    wtarget=w.*exp(2*t*real(z));wtarget=wtarget/sum(wtarget);
    verifyLessThan(tc,norm(S.cyclic.weights-wtarget),1e-12);
    verifyLessThan(tc,abs(S.finite.moments(2)-trace(H^2)/n),4e-10);
    V=schur_hessenberg_flow(n,t,z,w,'cmv','inverse');
    verifyLessThan(tc,norm(V'*V-eye(n),'fro'),2e-10);
    verifyLessThan(tc,norm(sort(angle(eig(V)))- ...
        sort(angle(eig(H)))),5e-8);
    F=schur_hessenberg_flow(n,t,z,w,'hessenberg','qr','operator');
    X=randn(n,2)+1i*randn(n,2);
    verifyLessThan(tc,norm(F(X,'notransp')-H*X,'fro'),2e-8);
end
end

function test_schur_3cycle_exact_weights_real_subclass(tc)
z=[1;exp(2*pi*1i/3);exp(-2*pi*1i/3)];
w=ones(3,1)/3;t=.7;
S=schur_hessenberg_flow(3,t,z,w,'hessenberg','qr','spectral');
target=[exp(3*t);1;1]/(exp(3*t)+2);
verifyLessThan(tc,norm(S.cyclic.weights-target),2e-13);
H=schur_hessenberg_flow(3,t,z,w);
verifyLessThan(tc,norm(imag(H),'fro'),2e-10);
verifyLessThan(tc,norm(H'*H-eye(3),'fro'),2e-10);
P=schur_hessenberg_flow(3,0,z,w,'hessenberg','qr','parameters');
verifyLessThan(tc,norm(P.hessenberg-P.initial_hessenberg,'fro'),1e-12);
end
