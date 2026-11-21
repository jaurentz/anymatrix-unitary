function tests=test_electric_quantum_walk
%TEST_ELECTRIC_QUANTUM_WALK Finite electric walk, cyclic CMV and spectra.
tests=functiontests(localfunctions);
end

function test_independent_shift_coin_electric_field_and_unitarity(tc)
for n=[3 4 6 9]
    for omega=[0 1/4 (sqrt(5)-1)/2]
        for theta=[0 .41 pi/2]
            phase=.13;beta=.27;
            U=electric_quantum_walk(n,omega,theta,phase,beta);
            Ref=direct_shift_coin_walk(n,omega,theta,phase,beta);
            D=electric_quantum_walk(n,omega,theta,phase,beta,'model');
            F=electric_quantum_walk(n,omega,theta,phase,beta,'factors');
            verifyEqual(tc,size(U),[2*n 2*n]);
            verifyTrue(tc,issparse(U));
            verifyLessThan(tc,norm(U'*U-speye(2*n),'fro'), ...
                1e-11*n);
            verifyLessThan(tc,norm(full(U-Ref),'fro'),1e-11*n);
            verifyLessThan(tc,norm(full(D.U-Ref),'fro'),1e-11*n);
            verifyLessThan(tc,norm(full(D.D*D.coin*D.S-U),'fro'), ...
                1e-11*n);
            verifyLessThan(tc,norm(full(F.D*F.L_even*F.M_odd-U),'fro'), ...
                1e-11*n);
            verifyLessThan(tc,norm(full(F.L_even'*F.L_even-speye(2*n)),'fro'), ...
                1e-11*n);
            verifyLessThan(tc,norm(full(F.M_odd'*F.M_odd-speye(2*n)),'fro'), ...
                1e-11*n);
            verifyLessThan(tc,norm(full(F.M_odd-F.M_odd'),'fro'), ...
                1e-11*n);
            verifyEqual(tc,F.standard_cyclic_cmv,false);
            verifyLessThan(tc,norm(abs(diag(F.D))-ones(2*n,1)),1e-13);
        end
    end
end
end

function test_exact_finite_ring_cmv_gauge_and_monodromy(tc)
% Independent original W from direct shift/coin construction was tested
% above; here the separately built cyclic CMV model must be gauge similar.
for n=[3 4 5 8 11]
    for omega=[0 1/4 (sqrt(5)-1)/2]
        for beta=[0 .37]
            for theta=[0 .43 pi/2]
                phase=.17;
                U=electric_quantum_walk(n,omega,theta,phase,beta);
                G=electric_quantum_walk(n,omega,theta,phase,beta, ...
                    'cmv_gauge');
                verifyTrue(tc,issparse(G.Lambda));
                verifyTrue(tc,issparse(G.L_even));
                verifyTrue(tc,issparse(G.M_odd));
                verifyTrue(tc,issparse(G.U_cmv));
                verifyLessThan(tc,norm(full(G.Lambda'*G.Lambda ...
                    -speye(2*n)),'fro'),5e-12*n);
                verifyLessThan(tc,norm(full(G.L_even'*G.L_even ...
                    -speye(2*n)),'fro'),5e-12*n);
                verifyLessThan(tc,norm(full(G.Lambda'*U*G.Lambda ...
                    -G.U_cmv),'fro'),2e-11*n);
                verifyLessThan(tc,norm(full(G.U_cmv ...
                    -G.global_phase*G.L_even*G.M_odd),'fro'), ...
                    2e-11*n);
                verifyLessThan(tc,abs(G.global_phase^n-G.monodromy), ...
                    2e-11*n);
                verifyLessThan(tc,norm(abs(G.alpha_even) ...
                    -sin(theta)*ones(n,1)),2e-11*n);
                verifyEqual(tc,G.alpha_odd,zeros(n,1));
                verifyLessThan(tc,norm(full(G.M_odd ...
                    -independent_odd_swap(n,beta)),'fro'),1e-12*n);
                lam=eig(full(U));mu=eig(full(G.U_cmv));
                verifyLessThan(tc,max(arrayfun(@(z)min(abs(mu-z)), ...
                    lam)),2e-9);
                for k=0:4
                    verifyLessThan(tc,abs(trace(full(U)^k) ...
                       -trace(full(G.U_cmv)^k)),2e-9*n);
                end
            end
        end
    end
end
% A nonzero electric field with unit monodromy admits an unscaled
% same-twist cyclic CMV; generic monodromy needs the scalar g.
n=5;omega=1/4;phi=1/2;
G=electric_quantum_walk(n,omega,.48,phi,.31,'cmv_gauge');
verifyLessThan(tc,abs(G.global_phase-1),1e-12);
verifyLessThan(tc,abs(G.monodromy-1),1e-12);
verifyLessThan(tc,norm(full(G.U_cmv-G.L_even*G.M_odd), ...
    'fro'),1e-11*n);
G=electric_quantum_walk(6,.271,.48,.17,.31,'cmv_gauge');
verifyGreaterThan(tc,abs(G.monodromy-1),.1);
verifyGreaterThan(tc,norm(full(G.U_cmv-G.L_even*G.M_odd), ...
    'fro'),.1);
% The sparse gauge and transformed CMV act on large multi-RHS arrays,
% without forming the default explicit walk (which is size-limited).
n=201;omega=(sqrt(5)-1)/2;
G=electric_quantum_walk(n,omega,.41,.09,-.37,'cmv_gauge');
O=electric_quantum_walk(n,omega,.41,.09,-.37,'operator');
X=reshape(sin((1:6*n)*.19),2*n,3) ...
    +1i*reshape(cos((1:6*n)*.21),2*n,3);
Y=G.Lambda*(G.U_cmv*(G.Lambda'*X));
verifyLessThan(tc,norm(O(X,'notransp')-Y,'fro'), ...
    5e-10*norm(X,'fro'));
verifyLessThan(tc,norm(O(X,'ctransp') ...
    -G.Lambda*(G.U_cmv'*(G.Lambda'*X)),'fro'), ...
    5e-10*norm(X,'fro'));
end

function test_standard_zero_field_cyclic_CMV_and_formal_skew_shift(tc)
for n=[3 5 8]
    for theta=[.21 .7 pi/2]
        E=electric_quantum_walk(n,0,theta,0,0,'factors');
        U=electric_quantum_walk(n,0,theta,0,0);
        verifyTrue(tc,E.standard_cyclic_cmv);
        b=sin(theta);a=cos(theta);
        if theta==pi/2,a=0;end
        thetaEven=[b,a;a,-b];
        expectedEven=kron(speye(n),sparse(thetaEven));
        expectedOdd=independent_odd_swap(n,0);
        verifyLessThan(tc,norm(full(E.L_even-expectedEven),'fro'), ...
            1e-13*n);
        verifyLessThan(tc,norm(full(E.M_odd-expectedOdd),'fro'), ...
            1e-13*n);
        verifyLessThan(tc,norm(full(E.D-speye(2*n)),'fro'), ...
            1e-13*n);
        verifyLessThan(tc,norm(full(U-expectedEven*expectedOdd),'fro'), ...
            2e-12*n);
        verifyLessThan(tc,norm(E.formal_alpha_even-b),2e-12*n);
        verifyEqual(tc,E.formal_alpha_odd,zeros(n,1));
        lambdas=eig(full(U));
        reference=eig(full(expectedEven*expectedOdd));
        verifyLessThan(tc,max(arrayfun(@(z)min(abs(reference-z)), ...
            lambdas)),2e-9);
        for k=[1 2 3]
            verifyLessThan(tc,abs(trace(full(U)^k) ...
              -trace(full(expectedEven*expectedOdd)^k)),1e-9);
        end
    end
end
for n=[5 9]
    omega=(sqrt(5)-1)/2;phi=.17;
    E=electric_quantum_walk(n,omega,.51,phi,.25,'factors');
    j=(0:n-1)';
    exact=sin(.51)*exp(2i*pi*mod(j.^2*omega+2*j*phi,1));
    verifyLessThan(tc,norm(E.formal_alpha_even-exact),1e-13);
    verifyEqual(tc,E.formal_alpha_odd,zeros(n,1));
    verifyFalse(tc,E.standard_cyclic_cmv);
end
end

function test_exact_transparent_coin_two_cycles_and_moments(tc)
for n=[3 4 7]
    for omega=[0 .25 (sqrt(5)-1)/2]
        beta=.37;
        U=electric_quantum_walk(n,omega,0,.14,beta);
        S=electric_quantum_walk(n,omega,0,.14,beta,'spectral');
        verifyEqual(tc,S.finite.kind,'exact');
        z=S.finite.eigenvalues();
        lam=eig(full(U));
        verifyEqual(tc,numel(z),2*n);
        verifyLessThan(tc,max(arrayfun(@(x)min(abs(lam-x)),z)), ...
            4e-9);
        for k=[-4 -1 0 1 2 3 4]
            verifyLessThan(tc,abs(S.finite.moments(k)-mean(z.^k)), ...
                2e-12);
            verifyLessThan(tc,abs(trace(full(U)^abs(k))/(2*n) ...
              -S.finite.moments(abs(k))),2e-9);
        end
        Q=U^n;
        pm=exp(1i*(2*pi*(n*.14+omega*n*(n-1)/2)-beta));
        pp=exp(1i*(2*pi*(n*.14+omega*n*(n-1)/2)+beta));
        verifyLessThan(tc,norm(full(Q(1:2:end,1:2:end) ...
             -pm*speye(n)),'fro'),2e-10*n);
        verifyLessThan(tc,norm(full(Q(2:2:end,2:2:end) ...
             -pp*speye(n)),'fro'),2e-10*n);
    end
end
end

function test_reflecting_bond_localization_exact_limit(tc)
for n=[5 7 10]
    omega=(sqrt(5)-1)/2;
    U=electric_quantum_walk(n,omega,pi/2,.17,.31);
    S=electric_quantum_walk(n,omega,pi/2,.17,.31,'spectral');
    verifyEqual(tc,S.finite.kind,'exact');
    z=S.finite.eigenvalues();
    [V,lam]=eig(full(U));
    verifyLessThan(tc,max(arrayfun(@(x)min(abs(diag(lam)-x)),z)), ...
        5e-9);
    % Each reflecting bond is an independent 2-by-2 unitary.
    % Every eigenvector has support on one bond: IPR = 1/2,
    % except for accidental degeneracies (not present here).
    for j=1:2*n
        v=V(:,j);
        ipr=sum(abs(v).^4)/sum(abs(v).^2)^2;
        verifyGreaterThan(tc,ipr,.49);
    end
    for k=[0 1 2 3 5]
        verifyLessThan(tc,abs(trace(full(U)^k)/(2*n) ...
          -S.finite.moments(k)),2e-10);
    end
end
end

function test_matrix_free_forward_adjoint_transpose_multirhs(tc)
for n=[4 9 30 72]
    omega=(sqrt(5)-1)/2;
    O=electric_quantum_walk(n,omega,.47,.08,.31,'operator');
    verifyEqual(tc,O([],'size'),[2*n 2*n]);
    X=reshape(sin(.41*(1:6*n)),2*n,3) ...
      +1i*reshape(cos(.25*(1:6*n)),2*n,3);
    Y=O(X,'notransp');
    verifyLessThan(tc,norm(O(Y,'ctransp')-X,'fro'),2e-12*n);
    verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')),2e-12*n);
    if n<=30
        U=electric_quantum_walk(n,omega,.47,.08,.31);
        verifyLessThan(tc,norm(Y-U*X,'fro'),2e-12*n);
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            2e-12*n);
        verifyLessThan(tc,norm(O(X,'transp')-U.'*X,'fro'), ...
            2e-12*n);
        verifyGreaterThan(tc,norm(O(X,'transp') ...
              -O(X,'ctransp'),'fro'),.1);
    end
    verifyError(tc,@()O(zeros(2*n-1,1),'notransp'), ...
        'um_electric_walk_apply:Dimension');
    verifyError(tc,@()O(X,'oops'),'um_electric_walk_apply:Flag');
end
end

function test_rational_irrational_phase_windows_and_global_shift(tc)
n=12;theta=.33;
R=electric_quantum_walk(n,1/4,theta,0,0,'factors');
d=full(diag(R.D));
verifyLessThan(tc,norm(d(1:8)-d(9:16)),1e-11);
A=electric_quantum_walk(n,(sqrt(5)-1)/2,theta,0,0,'factors');
z=full(diag(A.D));omega=(sqrt(5)-1)/2;
verifyGreaterThan(tc,abs(z(end)-z(1)),.1);
U=electric_quantum_walk(n,omega,theta,.08,.12);
V=electric_quantum_walk(n,omega,theta,.08+.21,.12);
verifyLessThan(tc,norm(full(V-exp(2i*pi*.21)*U),'fro'), ...
    2e-11*n);
Unknown=electric_quantum_walk(n,omega,theta,.08,.12,'spectral');
verifyEqual(tc,Unknown.finite.kind,'unknown');
verifyEmpty(tc,Unknown.finite.eigenvalues);
verifyEmpty(tc,Unknown.finite.moments);
end

function test_inventory_validation_and_large_operator(tc)
[W,p]=electric_quantum_walk();
verifyEmpty(tc,W);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyTrue(tc,any(strcmp(p,'sparse')));
verifyEqual(tc,size(electric_quantum_walk(4)),[8 8]);
verifyError(tc,@()electric_quantum_walk(2), ...
    'electric_quantum_walk:Order');
verifyError(tc,@()electric_quantum_walk(3,[],pi), ...
    'electric_quantum_walk:Parameters');
verifyError(tc,@()electric_quantum_walk(4,[],[],[],[],'bad'), ...
    'electric_quantum_walk:Mode');
verifyError(tc,@()electric_quantum_walk(2049), ...
    'electric_quantum_walk:ExplicitOrder');
O=electric_quantum_walk(2050,[],[],[],[],'operator');
verifyEqual(tc,O([],'size'),[4100 4100]);
x=complex(cos((1:4100)'*.03),sin((1:4100)'*.07));
verifyLessThan(tc,norm(O(O(x,'notransp'),'ctransp')-x), ...
    1e-9*norm(x));
end

function U=direct_shift_coin_walk(n,omega,theta,phase,beta)
% Independent physical Q*S realization: site ordered (-,+).
a=cos(theta);b=sin(theta);
if theta==pi/2,a=0;end
C=kron(speye(n),sparse([a,b;-b,a]));
S=sparse(2*n,2*n);
for j=0:n-1
    colm=2*j+1;colmto=2*mod(j-1,n)+1;
    colp=2*j+2;colpto=2*mod(j+1,n)+2;
    S(colmto,colm)=exp(-1i*beta*(j==0));
    S(colpto,colp)=exp(1i*beta*(j==n-1));
end
phases=exp(2i*pi*mod(phase+omega*(0:n-1)',1));
D=spdiags(repelem(phases,2),0,2*n,2*n);
U=D*C*S;
end

function M=independent_odd_swap(n,beta)
M=sparse(2*n,2*n);
for j=0:n-1
    plus=2*j+2;
    nextminus=2*mod(j+1,n)+1;
    phi=exp(1i*beta*(j==n-1));
    M(nextminus,plus)=phi;
    M(plus,nextminus)=conj(phi);
end
end
