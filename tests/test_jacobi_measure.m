function tests = test_jacobi_measure
%TEST_JACOBI_MEASURE Independent Jacobi, Szego and Geronimus tests.
tests = functiontests(localfunctions);
end

function test_normalization(testCase)
for ab = [0 0; -0.5 -0.5; 0.5 0.5; 0.3 1.6; 2.4 0.7]'
    a=ab(1); b=ab(2);
    mu=szego_jacobi_measure(a,b);
    mass = integral(mu.circle_density,-pi,pi,'AbsTol',5e-11,'RelTol',5e-11);
    verifyLessThan(testCase,abs(mass-1),2e-8);
    massReal=integral(mu.real_density,-2,2,'AbsTol',5e-11,'RelTol',5e-11);
    verifyLessThan(testCase,abs(massReal-1),2e-8);
    verifyEqual(testCase,mu.family,'jacobi_szego');
end
end

function test_explicit_and_inverse_geronimus(testCase)
parameters = [0 0; -0.5 -0.5; 0.5 0.5; 0.3 1.6; -0.6 0.4; 2.4 -0.7];
for j=1:size(parameters,1)
    a=parameters(j,1); b=parameters(j,2);
    mu=szego_jacobi_measure(a,b);
    alpha=szego_verblunsky(mu,24);
    verifyEqual(testCase,size(alpha),[24,1]);
    verifyTrue(testCase,all(abs(alpha)<1));
    [d,s]=js_jacobi_recurrence(a,b,13);
    reconstructed=js_inverse_geronimus(d,s,24);
    verifyLessThan(testCase,norm(alpha-reconstructed,inf),1e-9);
    [d2,s2]=js_forward_geronimus(alpha,12);
    verifyLessThan(testCase,norm(d2-d(1:12),inf),2e-11);
    verifyLessThan(testCase,norm(s2-s(1:12),inf),2e-11);
end
end

function test_classical_specializations(testCase)
verifyLessThan(testCase,norm(szego_verblunsky(szego_jacobi_measure(-0.5,-0.5),21)),1e-14);
leg=szego_verblunsky(szego_jacobi_measure(0,0),18);
verifyLessThan(testCase,norm(leg(1:2:end)),1e-14);
verifyLessThan(testCase,norm(leg(2:2:end) + 1./(3:2:19)'),1e-13);
end

function test_exact_low_order_moments(testCase)
mu=szego_jacobi_measure(0,0);
% Legendre measure: E(cos(theta))=0 and E(cos(2theta))=-1/3.
m=js_jacobi_cosine_moments(mu.a,mu.b,4);
verifyEqual(testCase,m(1),1,'AbsTol',1e-13);
verifyEqual(testCase,m(2),0,'AbsTol',1e-13);
verifyEqual(testCase,m(3),-1/3,'AbsTol',1e-13);
end

function test_invalid_measure_parameters(testCase)
verifyTrue(testCase,raises(@()szego_jacobi_measure(-1,0)));
verifyTrue(testCase,raises(@()szego_jacobi_measure(0,-1)));
end

function tf=raises(f)
tf=false;
try
    f();
catch
    tf=true;
end
end
