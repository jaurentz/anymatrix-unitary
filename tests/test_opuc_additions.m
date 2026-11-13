function tests=test_opuc_additions
%TEST_OPUC_ADDITIONS Measure-first Bernstein-Szego and circular-beta models.
tests=functiontests(localfunctions);
end

function test_bernstein_szego_measure(testCase)
for a=[0 .25+0.4i -.6 .8]
    mu=szego_bernstein_measure(a);
    verifyEqual(testCase,mu.family,'bernstein_szego');
    mass=integral(mu.circle_density,-pi,pi,'AbsTol',1e-11);
    verifyLessThan(testCase,abs(mass-1),1e-10);
    first=integral(@(theta)exp(1i*theta).*mu.circle_density(theta), ...
        -pi,pi,'AbsTol',1e-11);
    verifyLessThan(testCase,abs(first-conj(a)),1e-10);
    alpha=szego_verblunsky(mu,18);
    verifyLessThan(testCase,abs(alpha(1)-a),1e-14);
    verifyEqual(testCase,alpha(2:end),zeros(17,1));
end
end

function test_bernstein_szego_matrices(testCase)
for n=[1 2 3 7 12]
    a=.35+.2i;mu=szego_bernstein_measure(a);
    H=szego_matrix(mu,n,'hessenberg',exp(1i*.35));
    for representation={'cmv','hessenberg','inverse_hessenberg','twisted'}
        U=bernstein_szego(n,a,representation{1},exp(1i*.35));
        verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),2e-12*n);
        for k=1:n
            verifyLessThan(testCase,abs(trace(U^k)-trace(H^k)),1e-9*n);
        end
        [z,w]=szego_quadrature(U);
        if n>1
            verifyLessThan(testCase,abs(sum(w.*z)-conj(a)),1e-10);
        end
    end
end
end

function test_circular_beta_seed_and_structure(testCase)
s0=rng;
for n=[1 2 3 8 19]
    for beta=[.5 1 2 4 8]
        U=circular_beta(n,beta,111);
        V=circular_beta(n,beta,111);
        verifyEqual(testCase,U,V);
        verifyTrue(testCase,issparse(U));
        verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),2e-11);
        verifyLessThanOrEqual(testCase,nnz(U),5*n);
    end
end
verifyEqual(testCase,rng,s0);
end

function test_circular_beta_radial_law(testCase)
% Elementary numerical check of the coefficient law, not a test of
% circular-ensemble level spacing in a finite sample.
rngState=rng; restore=onCleanup(@()rng(rngState)); %#ok<NASGU>
rng(12,'twister');
N=1800;n=6;beta=2;
r=zeros(N,1);
for k=1:N
    alpha=js_circular_beta_coeffs(n,beta);
    r(k)=abs(alpha(2))^2; % alpha_1, Beta(1,beta*(n-2)/2)
end
expected=1/(1+beta*(n-2)/2);
verifyLessThan(testCase,abs(mean(r)-expected),.025);
end

function test_bad_opuc_parameters(testCase)
verifyError(testCase,@()szego_bernstein_measure(1),'js_bernstein_measure:Alpha');
verifyTrue(testCase,throws_error(@()circular_beta(4,0)));
end

function yes=throws_error(f)
yes=false;
try
    f();
catch
    yes=true;
end
end
