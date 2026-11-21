function tests=test_fourier_laguerre_ball
tests=functiontests(localfunctions);
end

function test_square_radial_zonal_orthogonality_and_operator(tc)
for P=[1 2 4 8]
    for L=[1 2 4 6]
        n=P*L;
        U=fourier_laguerre_ball(P,L,.8);
        verifyEqual(tc,size(U),[n n]);
        verifyTrue(tc,isreal(U));
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),4e-11*n);
        par=fourier_laguerre_ball(P,L,.8,'parameters');
        verifyLessThan(tc,norm(U-kron(par.angular_transform, ...
            par.radial_transform),'fro'),2e-12*n);
        verifyLessThan(tc,norm(par.angular_transform'* ...
            par.angular_transform-eye(L),'fro'),1e-11*n);
        verifyLessThan(tc,norm(par.radial_transform'* ...
            par.radial_transform-eye(P),'fro'),1e-11*n);
        op=fourier_laguerre_ball(P,L,.8,'operator');
        X=randn(n,3)+1i*randn(n,3);
        verifyEqual(tc,op([],'size'),[n n]);
        verifyLessThan(tc,norm(op(X,'notransp')-U*X,'fro'),4e-11*n);
        verifyLessThan(tc,norm(op(X,'transp')-U'*X,'fro'),4e-11*n);
        verifyLessThan(tc,norm(op(op(X,'notransp'),'ctransp')-X,'fro'), ...
            4e-11*n);
    end
end
end

function test_independent_radial_laguerre_evaluation(tc)
for P=[1 2 3 6 10]
    tau=1.7;
    par=fourier_laguerre_ball(P,3,tau,'parameters');
    x=par.laguerre_nodes;
    w=par.radial_quadrature_weights;
    Q=zeros(P);
    old=ones(P,1);prev=zeros(P,1);
    for p=0:P-1
        if p==0
            poly=old;
        elseif p==1
            poly=3-x;
        else
            poly=((2*p+1-x).*old-(p+1)*prev)/p;
            % Recurrence: L_p^2=((2(p-1)+3-x)L_(p-1)^2
            %                  -((p-1)+2)L_(p-2)^2)/p.
        end
        Q(:,p+1)=sqrt(w).*sqrt(exp(gammaln(p+1)-gammaln(p+3))) ...
            .*tau^(-3/2).*exp(-x/2).*poly;
        prev=old;old=poly;
    end
    verifyLessThan(tc,norm(Q-par.radial_transform,'fro'),1e-8*P);
    verifyLessThan(tc,abs(sum(par.laguerre_weighted_weights)-2),1e-12);
    verifyGreaterThan(tc,min(par.radial_nodes),0);
end
end

function test_independent_angular_quadrature(tc)
for L=[1 2 3 7]
    par=fourier_laguerre_ball(3,L,1,'parameters');
    x=par.cos_theta_nodes;w=par.theta_weights;
    A=zeros(L);
    for ell=0:L-1
        Pl=legendre(ell,x.');
        A(:,ell+1)=sqrt(w*(2*ell+1)/2).*Pl(1,:).';
    end
    verifyLessThan(tc,norm(A-par.angular_transform,'fro'),3e-11*L);
    verifyLessThan(tc,abs(sum(w)-2),3e-13);
end
end

function test_full_spherical_sampling_rectangular_isometry(tc)
for P=[1 2 3]
    for L=[1 2 4]
        A=fourier_laguerre_ball(P,L,1,'sampling');
        verifyEqual(tc,size(A),[P*L*(2*L-1),P*L^2]);
        verifyLessThan(tc,norm(A'*A-eye(P*L^2),'fro'), ...
            4e-10*P*L^2);
        if L>1
            verifyGreaterThan(tc,size(A,1),size(A,2));
        end
    end
end
end

function test_physical_tau_and_unknown_spectrum(tc)
P=5;L=4;
A=fourier_laguerre_ball(P,L,.5);
B=fourier_laguerre_ball(P,L,3.1);
verifyEqual(tc,A,B,'AbsTol',3e-13);
a=fourier_laguerre_ball(P,L,.5,'parameters');
b=fourier_laguerre_ball(P,L,3.1,'parameters');
verifyLessThan(tc,norm(b.radial_nodes-(3.1/.5)*a.radial_nodes),3e-12);
S=fourier_laguerre_ball(P,L,.5,'spectral');
verifyEqual(tc,S.family,'fourier_laguerre_ball');
verifyEqual(tc,S.dimension,P*L);
verifyEqual(tc,S.finite.kind,'unknown');
verifyError(tc,@()fourier_laguerre_ball(P,L,1,'invalid'), ...
    'fourier_laguerre_ball:Mode');
verifyError(tc,@()fourier_laguerre_ball(30,30,1,'sampling'), ...
    'fourier_laguerre_ball:SamplingSize');
end
