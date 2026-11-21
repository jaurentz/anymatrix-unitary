function tests=test_discrete_fractional_fourier
%TEST_DISCRETE_FRACTIONAL_FOURIER Exact Candan parity-ordered FrFT.
tests=functiontests(localfunctions);
end

function test_dft_endpoint_and_even_index_gap(tc)
for n=[1 2 3 4 5 6 7 8 9 12 15 16 20 23 24]
    F=discrete_fourier(n);
    B=discrete_fractional_fourier(n,1,'basis');
    V=B.V;k=B.hermite_indices;
    if mod(n,2)==0
        verifyEqual(tc,k(:),[(0:n-2)';n]);
    else
        verifyEqual(tc,k(:),(0:n-1)');
    end
    verifyEqual(tc,size(V),[n n]);
    verifyLessThan(tc,norm(imag(V),'fro'),4e-14*n);
    verifyLessThan(tc,norm(V'*V-eye(n),'fro'),3e-11*n);
    verifyLessThan(tc,norm(F*V-V*diag(exp(-.5i*pi*k)),'fro'), ...
        5e-10*n);
    verifyLessThan(tc,norm(discrete_fractional_fourier(n,1)-F,'fro'), ...
        6e-10*n);
    verifyLessThan(tc,norm(discrete_fractional_fourier(n,0)-eye(n),'fro'), ...
        5e-11*n);
    verifyLessThan(tc,norm(discrete_fractional_fourier(n,2)-F^2,'fro'), ...
        6e-10*n);
    verifyLessThan(tc,norm(discrete_fractional_fourier(n,4)-eye(n),'fro'), ...
        6e-10*n);
end
end

function test_commuting_oscillator_and_degeneracy_parity(tc)
for n=[1 2 4 8 12 16 20 25]
    F=discrete_fourier(n);
    B=discrete_fractional_fourier(n,.2,'basis');
    S=discrete_fractional_fourier(n,.4,'oscillator');
    verifyLessThan(tc,norm(S-S','fro'),2e-13*n);
    verifyLessThan(tc,norm(S*F-F*S,'fro'),1e-11*n);
    verifyLessThan(tc,norm(B.V'*S*B.V-diag(B.oscillator_eigenvalues),'fro'), ...
        3e-10*n);
    for j=1:n
        i=(0:n-1)';
        paritySample=B.V(mod(-i,n)+1,j);
        verifyLessThan(tc,norm(paritySample-B.parity(j)*B.V(:,j)), ...
            2e-12*n);
    end
    for j=1:n
        [~,pivot]=max(abs(B.V(:,j)));
        verifyGreaterThanOrEqual(tc,B.V(pivot,j),-2e-13);
    end
    if mod(n,4)==0
        ev=B.oscillator_eigenvalues;
        verifyEqual(tc,sum(abs(ev+4)<1e-9),2);
        zeroidx=find(abs(ev+4)<1e-9);
        verifyEqual(tc,sort(B.parity(zeroidx)),[-1;1]);
    end
end
end

function test_additivity_adjoint_periodicity_and_hamiltonian(tc)
rng(47);
for n=[1 3 4 7 8 11 16]
    for a=[-.2 .3 .87 1.42]
        b=.23;
        A=discrete_fractional_fourier(n,a);
        B=discrete_fractional_fourier(n,b);
        AB=discrete_fractional_fourier(n,a+b);
        invA=discrete_fractional_fourier(n,-a);
        four=discrete_fractional_fourier(n,a+4);
        verifyLessThan(tc,norm(A'*A-eye(n),'fro'),2e-10*n);
        verifyLessThan(tc,norm(A*B-AB,'fro'),4e-10*n);
        verifyLessThan(tc,norm(A'-invA,'fro'),2e-10*n);
        verifyLessThan(tc,norm(A-four,'fro'),4e-10*n);
        H=discrete_fractional_fourier(n,a,'hamiltonian');
        verifyLessThan(tc,norm(H-H','fro'),2e-12*n);
        verifyLessThan(tc,norm(expm(-.5i*pi*a*H)-A,'fro'), ...
            2e-9*n);
    end
end
end

function test_exact_trace_spectrum_and_moments(tc)
for n=[1 2 4 5 8 13 16]
    for a=[.2 1 2]
        Q=discrete_fractional_fourier(n,a);
        S=discrete_fractional_fourier(n,a,'spectral');
        verifyEqual(tc,S.family,'discrete_fractional_fourier');
        verifyEqual(tc,S.finite.kind,'exact');
        kk=S.hermite_indices;
        lambda=exp(-.5i*pi*a*kk);
        verifyLessThan(tc,norm(S.finite.eigenvalues()-lambda),2e-14*n);
        verifyLessThan(tc,abs(sum(S.finite.atoms.weights)-1),5e-15);
        for p=-4:4
            verifyLessThan(tc,abs(S.finite.moments(p)- ...
                mean(lambda.^p)),3e-14*n);
            verifyLessThan(tc,abs(S.finite.moments(p)- ...
                trace(Q^p)/n),2e-9*n);
        end
    end
end
end

function test_reproducibility_and_domain(tc)
A=discrete_fractional_fourier(16,.373,'basis');
B=discrete_fractional_fourier(16,.373,'basis');
verifyEqual(tc,A.V,B.V,'AbsTol',0);
verifyError(tc,@()discrete_fractional_fourier(7,.3,'operator'), ...
    'discrete_fractional_fourier:Mode');
verifyError(tc,@()discrete_fractional_fourier(4,.3,'bad'), ...
    'discrete_fractional_fourier:Mode');
verifyError(tc,@()discrete_fractional_fourier(257), ...
    'MATLAB:notLessEqual');
end
