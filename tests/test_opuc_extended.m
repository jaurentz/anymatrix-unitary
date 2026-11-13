function tests=test_opuc_extended
%TEST_OPUC_EXTENDED Real/complex OPUC measures and finite unitary matrices.
tests=functiontests(localfunctions);
end

function test_explicit_verblunsky_sequences(testCase)
a=.4;eta=.7;b=a+1i*eta;
mu=szego_circular_jacobi_measure(a,eta);
alpha=szego_verblunsky(mu,12);
expected=zeros(12,1);
r=1;
for k=0:11
    r=r*(b+k)/(1+conj(b)+k);
    expected(k+1)=-r;
end
verifyLessThan(testCase,max(abs(alpha-expected)),2e-14);
verifyLessThan(testCase,max(abs(szego_verblunsky( ...
    szego_circular_jacobi_measure(a,0),12) + ...
    a./(a+(1:12)'))),2e-14);

q=.64;
verifyLessThan(testCase,max(abs(szego_verblunsky( ...
    szego_rogers_measure(q),5) - ...
    [sqrt(q);-q;q^(3/2);-q^2;q^(5/2)])),1e-14);

verifyEqual(testCase,szego_verblunsky( ...
    szego_geronimus_measure(.3+.2i),4),repmat(.3+.2i,4,1));
verifyEqual(testCase,szego_verblunsky( ...
    szego_periodic_measure([.2;.3i;-.4]),8), ...
    [.2;.3i;-.4;.2;.3i;-.4;.2;.3i]);
verifyEqual(testCase,szego_verblunsky( ...
    szego_bernstein_finite_measure([.2;.3i]),5),[.2;.3i;0;0;0]);
end

function test_explicit_circle_weights(testCase)
q=.64;
mu=szego_rogers_measure(q);
mass=integral(mu.circle_density,-pi,pi,'AbsTol',1e-12);
m1=integral(@(t)exp(1i*t).*mu.circle_density(t),-pi,pi,'AbsTol',1e-12);
m2=integral(@(t)exp(2i*t).*mu.circle_density(t),-pi,pi,'AbsTol',1e-12);
verifyLessThan(testCase,abs(mass-1),1e-10);
verifyLessThan(testCase,abs(m1-sqrt(q)),1e-10);
verifyLessThan(testCase,abs(m2-q^2),1e-10);

a=.4;eta=.7;b=a+1i*eta;
mu=szego_circular_jacobi_measure(a,eta);
mass=integral(mu.circle_density,0,2*pi,'AbsTol',1e-11);
m1=integral(@(t)exp(1i*t).*mu.circle_density(t),0,2*pi,'AbsTol',1e-11);
verifyLessThan(testCase,abs(mass-1),1e-9);
verifyLessThan(testCase,abs(m1-conj(-b/(1+conj(b)))),1e-9);

alpha=[.25;.2i;-.17];
mu=szego_bernstein_finite_measure(alpha);
mass=integral(mu.circle_density,-pi,pi,'AbsTol',1e-11);
m1=integral(@(t)exp(1i*t).*mu.circle_density(t),-pi,pi,'AbsTol',1e-11);
verifyLessThan(testCase,abs(mass-1),1e-9);
verifyLessThan(testCase,abs(m1-conj(alpha(1))),1e-9);
% The general finite family reduces exactly to the existing Poisson weight.
a=.2+.3i;
old=szego_bernstein_measure(a);
new=szego_bernstein_finite_measure(a);
t=linspace(-pi,pi,63);
verifyLessThan(testCase,max(abs(old.circle_density(t)-new.circle_density(t))),1e-13);
end

function test_new_matrices_and_four_representations(testCase)
forms={'cmv','hessenberg','inverse_hessenberg','twisted'};
cases={szego_circular_jacobi_measure(.4,.7), ...
       szego_rogers_measure(.6), ...
       szego_geronimus_measure(.2+.15i), ...
       szego_periodic_measure([.2;.3i]), ...
       szego_bernstein_finite_measure([.25;.2i;-.1])};
for i=1:numel(cases)
    mu=cases{i};
    for n=[1 2 3 8 16]
        beta=exp(.31i);
        H=szego_matrix(mu,n,'hessenberg',beta);
        for j=1:numel(forms)
            U=szego_matrix(mu,n,forms{j},beta);
            verifyTrue(testCase,issparse(U));
            verifyEqual(testCase,size(U),[n,n]);
            verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),2e-12*n);
            verifyLessThan(testCase,abs(trace(U)-trace(H)),2e-10*n);
            [z,w]=szego_quadrature(U);
            verifyLessThan(testCase,abs(sum(w)-1),1e-11);
            if n>1
                a=szego_verblunsky(mu,1);
                verifyLessThan(testCase,abs(sum(w.*z)-conj(a)),2e-10);
            end
        end
    end
end
end

function test_public_generator_delegation(testCase)
n=11;beta=exp(.27i);
U={circular_jacobi(n,.3,.7,'cmv',beta), ...
   rogers_szego(n,.4,'cmv',beta), ...
   geronimus_cmv(n,.1+.2i,'cmv',beta), ...
   periodic_cmv(n,[.1;.2i],'cmv',beta), ...
   bernstein_szego_finite(n,[.1;.2i],'cmv',beta)};
for k=1:numel(U)
    verifyLessThan(testCase,norm(full(U{k}'*U{k}-speye(n)),'fro'),2e-11);
    verifyLessThanOrEqual(testCase,nnz(U{k}),5*n);
end
for generator={'circular_jacobi','rogers_szego','geronimus_cmv', ...
        'periodic_cmv','bernstein_szego_finite'}
    [empty,props]=feval(generator{1});
    verifyEmpty(testCase,empty);
    verifyTrue(testCase,any(strcmp(props,'unitary')));
end
% Real coefficients and a real closure imply a real orthogonal matrix.
for mu={szego_rogers_measure(.4), szego_geronimus_measure(.3), ...
        szego_periodic_measure([.2;-.1]), ...
        szego_bernstein_finite_measure([.2;.1])}
    R=szego_matrix(mu{1},12,'cmv',-1);
    verifyTrue(testCase,isreal(R));
    verifyLessThan(testCase,norm(full(R'*R-speye(12)),'fro'),1e-11);
end
end

function test_parameter_rejection(testCase)
verifyError(testCase,@()szego_rogers_measure(1),'MATLAB:expectedLess');
verifyError(testCase,@()szego_geronimus_measure(1), ...
    'szego_geronimus_measure:Alpha');
verifyError(testCase,@()szego_periodic_measure([.2;1]), ...
    'szego_periodic_measure:Alpha');
verifyError(testCase,@()szego_bernstein_finite_measure([.1;1]), ...
    'szego_bernstein_finite_measure:Alpha');
verifyError(testCase,@()circular_jacobi(5,.2,.1,'cmv',0.7), ...
    'js_core_product:Beta');
end
