function tests=test_spectral_examples
%TEST_SPECTRAL_EXAMPLES Execute existing examples and check band formulas.
tests=functiontests(localfunctions);
end

function test_quickstart_example(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
original=pwd;
cleanup=onCleanup(@()cd(original)); %#ok<NASGU>
cd(root);
run(fullfile(root,'examples','quickstart.m'));
verifyEqual(testCase,size(D),[96,96]);
verifyLessThan(testCase,norm(full(D'*D-speye(96)),'fro'),1e-11);
verifyEqual(testCase,Pop([],'size'),[48,48]);
verifyLessThan(testCase,abs(sum(weights)-1),1e-12);
end

function test_showcase_without_plot(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
original=pwd;
cleanup=onCleanup(@()cd(original)); %#ok<NASGU>
cd(root);
plotSpectralBands=false;
run(fullfile(root,'examples','showcase_new_families.m'));
verifyEqual(testCase,numel(phaseBands),56);
verifyEqual(testCase,size(D),[192,192]);
verifyEqual(testCase,size(P),[96,96]);
verifyEqual(testCase,size(S),[144,144]);
verifyTrue(testCase,any(inBandP) && any(~inBandP));
verifyGreaterThan(testCase,positiveEdgesP(1),0);
verifyLessThan(testCase,positiveEdgesP(2),pi);
verifyLessThanOrEqual(testCase,numBoundaryLevels,96);
verifyGreaterThan(testCase,rhoD(find(insideD,1)),0);
verifyTrue(testCase,all(rhoD(~insideD)==0));
verifyTrue(testCase,all(rhoT(~insideT)==0));
verifyTrue(testCase,islogical(plotSpectralBands) && ~plotSpectralBands);
end

function test_period_two_discriminant_from_transfer_matrices(testCase)
% Standard Szego transfer matrices, for real period-two coefficients.
a0=.48; a1=-.31;
rho0=sqrt(1-a0^2);rho1=sqrt(1-a1^2);
for theta=linspace(-pi,pi,19)
    z=exp(1i*theta);
    A0=[z,-a0;-a0*z,1]/rho0;
    A1=[z,-a1;-a1*z,1]/rho1;
    discr=trace(A1*A0)/z;
    formula=2*(cos(theta)+a0*a1)/(rho0*rho1);
    verifyLessThan(testCase,abs(discr-formula),2e-12);
end
% Constant real coefficients are a special case of period two.
a=.36;
edge=2*asin(a);
verifyLessThan(testCase,abs(cos(edge)-(1-2*a^2)),2e-14);
end

function test_dirac_density_and_cycle_limit(testCase)
m=.45;nu=sqrt(1-m^2);
edge=asin(m);
theta=[0;edge/2;pi/2;pi-edge/2;pi];
support=abs(cos(theta))<nu;
verifyEqual(testCase,support,[false;false;true;false;false]);
rho=abs(sin(theta(support)))./ ...
    (2*pi*sqrt(nu^2-cos(theta(support)).^2));
verifyLessThan(testCase,abs(rho-1/(2*pi*nu)),1e-14);
% At fixed t without circular wrapping, the normalized trace converges
% exponentially in n for this analytic periodic integrand.
n=96;t=.65;
T=cycle_propagator(n,t);
verifyLessThan(testCase,abs(trace(T)/n-besselj(0,2*t)),2e-12);
end
