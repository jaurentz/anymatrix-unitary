function tests=test_hermite_oscillator
%TEST_HERMITE_OSCILLATOR Exact Hermite-quadrature spectral unitary flow.
tests=functiontests(localfunctions);
end

function test_orthogonal_Gauss_Hermite_basis_and_quadrature(tc)
for n=[1 2 3 5 9 20]
    B=hermite_oscillator(n,.7,0,'basis');
    Q=B.Q;x=B.nodes;w=B.weights;
    verifyEqual(tc,size(Q),[n n]);
    verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),3e-12*n);
    verifyLessThan(tc,norm(Q-Q*0-gauss_hermite_transform(n),'fro'), ...
        1e-13);
    verifyGreaterThan(tc,min(w),0);
    verifyLessThan(tc,abs(sum(w)-1),3e-12);
    verifyEqual(tc,numel(x),n);
    verifyTrue(tc,all(diff(x)>0));
    verifyLessThan(tc,abs(sum(w.*x)),2e-12);
    if n>=2
        verifyLessThan(tc,abs(sum(w.*x.^2)-.5),2e-11);
        verifyLessThan(tc,norm(Q(:,2)-sqrt(2)*x.*Q(:,1)), ...
            3e-11);
    end
    if n>=3
        verifyLessThan(tc,abs(sum(w.*x.^4)-.75),2e-11);
        for k=2:n-1
            next=sqrt(2/(k))*x.*Q(:,k) ...
                -sqrt((k-1)/k)*Q(:,k-1);
            verifyLessThan(tc,norm(Q(:,k+1)-next),3e-9);
        end
    end
end
end

function test_unitarity_group_additivity_periodicity(tc)
for n=[1 2 3 6 10]
    for offset=[0 .5]
        for theta=[0 .27 -.8 pi/2 pi]
            U=hermite_oscillator(n,theta,offset);
            V=hermite_oscillator(n,-theta,offset);
            verifyEqual(tc,size(U),[n n]);
            verifyLessThan(tc,norm(U'*U-eye(n),'fro'),3e-12*n);
            verifyLessThan(tc,norm(U*V-eye(n),'fro'),3e-12*n);
            verifyLessThan(tc,norm(U'-V,'fro'),3e-12*n);
            verifyLessThan(tc,norm(hermite_oscillator(n,theta+2*pi,offset) ...
                -(-1)^(2*offset)*U,'fro'),4e-12*n);
            if theta==0
                verifyLessThan(tc,norm(U-eye(n),'fro'),3e-12*n);
            end
        end
        U=hermite_oscillator(n,.43,offset);
        V=hermite_oscillator(n,-.92,offset);
        W=hermite_oscillator(n,-.49,offset);
        verifyLessThan(tc,norm(U*V-W,'fro'),5e-12*n);
        verifyLessThan(tc,norm(U*V-V*U,'fro'),5e-12*n);
    end
end
end

function test_exact_finite_spectrum_and_trace_moments(tc)
for n=[1 2 3 8 16]
    for offset=[0 .5]
        for theta=[0 .36 pi/2 pi 2*pi]
            U=hermite_oscillator(n,theta,offset);
            S=hermite_oscillator(n,theta,offset,'spectral');
            verifyEqual(tc,S.family,'hermite_oscillator');
            verifyEqual(tc,S.dimension,n);
            verifyEqual(tc,S.finite.kind,'exact');
            z=S.finite.eigenvalues();
            energy=(0:n-1)'+offset;
            verifyLessThan(tc,norm(z-exp(-1i*theta*energy)), ...
                4e-13*n);
            verifyEqual(tc,S.oscillator_energies,energy);
            lam=eig(U);
            verifyLessThan(tc,max(arrayfun(@(x)min(abs(lam-x)),z)), ...
                2e-9);
            for k=[-5 -2 -1 0 1 2 5]
                exact=mean(exp(-1i*k*theta*energy));
                verifyLessThan(tc,abs(S.finite.moments(k)-exact), ...
                    2e-12);
                verifyLessThan(tc,abs(trace(U^abs(k))/n ...
                  -S.finite.moments(abs(k))),3e-10);
            end
        end
    end
end
end

function test_independent_harmonic_oscillator_Hamiltonian(tc)
for n=[2 4 7 12]
    for offset=[0 .5]
        H=hermite_oscillator(n,.1,offset,'hamiltonian');
        verifyLessThan(tc,norm(H-H','fro'),2e-12*n);
        verifyTrue(tc,isreal(H));
        verifyLessThan(tc,norm(sort(eig(H)) ...
            -((0:n-1)'+offset)),5e-11*n);
        for theta=[-.43 .1 .92]
            U=hermite_oscillator(n,theta,offset);
            verifyLessThan(tc,norm(U-expm(-1i*theta*H),'fro'), ...
                5e-11*n);
            verifyLessThan(tc,norm(H*U-U*H,'fro'),5e-11*n);
            delta=1e-5;
            approx=(hermite_oscillator(n,theta+delta,offset) ...
                -hermite_oscillator(n,theta-delta,offset))/(2*delta);
            verifyLessThan(tc,norm(approx+1i*H*U,'fro'), ...
                2e-7*n*n);
        end
    end
end
end

function test_not_a_DFT_compatible_discrete_FrFT(tc)
% Distinguish the strict GH quadrature-sample spectral model from
% Candan-Kutay-Ozaktas' DFT-compatible discrete Hermite basis.
for n=[3 5 8 12]
    U=hermite_oscillator(n,pi/2,0);
    F=discrete_fourier(n);
    verifyGreaterThan(tc,norm(U-F,'fro'),.1);
    verifyGreaterThan(tc,norm(U-conj(F),'fro'),.1);
    verifyGreaterThan(tc,norm(U-diag(diag(U)),'fro'),.1);
end
end

function test_properties_and_parameter_validation(tc)
[U,p]=hermite_oscillator();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyTrue(tc,any(strcmp(p,'complex')));
verifyTrue(tc,any(strcmp(p,'eigenvalues')));
verifyEqual(tc,size(hermite_oscillator(5)),[5 5]);
verifyError(tc,@()hermite_oscillator(0), ...
    'hermite_oscillator:Order');
verifyError(tc,@()hermite_oscillator(3,1i), ...
    'hermite_oscillator:Theta');
verifyError(tc,@()hermite_oscillator(3,.2,.37), ...
    'hermite_oscillator:Offset');
verifyError(tc,@()hermite_oscillator(3,.2,0,'invalid'), ...
    'hermite_oscillator:Mode');
end
