function tests=test_geronimus_spectral_measure
%TEST_GERONIMUS_SPECTRAL_MEASURE Analytic infinite e1 versus finite CMV measures.
tests=functiontests(localfunctions);
end

function test_infinite_real_and_complex_continuous_atom_split(tc)
cases={.8,-.4,.2+.3i,-.3+.6i,-.6-.1i,0};
theta=linspace(-pi,pi,32769);
for j=1:numel(cases)
    a=cases{j};
    mu=szego_geronimus_measure(a);
    rho=mu.circle_density(theta);
    verifyGreaterThanOrEqual(tc,min(rho(:)),-1e-13);
    atom=mu.point_mass.weights;
    if isempty(atom),atom=0;end
    mass=trapz(theta,rho)+atom;
    verifyLessThan(tc,abs(mass-1),2e-5);
    phase=mu.point_mass.angles;
    if isempty(phase),phase=0;end
    moment=trapz(theta,exp(1i*theta).*rho)+atom*exp(1i*phase);
    verifyLessThan(tc,abs(moment-conj(a)),2e-5);
    verifyLessThan(tc,abs(mu.ac_mass+atom-1),3e-15);
    if a~=0
        verifyEqual(tc,mu.circle_density(0),0);
        verifyEqual(tc,mu.circle_density(2*pi),0);
    end
end
mu=szego_geronimus_measure(.8);
verifyLessThan(tc,abs(mu.point_mass.weights-8/9),3e-15);
verifyLessThan(tc,abs(mu.ac_mass-1/9),3e-15);
verifyLessThan(tc,abs(mu.point_mass.angles),3e-15);
verifyEmpty(tc,szego_geronimus_measure(-.4).point_mass.weights);
verifyLessThan(tc,abs(szego_geronimus_measure(0).circle_density(0)-1/(2*pi)),1e-15);
end

function test_poisson_is_exact_analytic_schur_smoothing(tc)
theta=linspace(-pi,pi,401);
for a=[0 .8 -.4 .2+.3i]
    mu=szego_geronimus_measure(a);
    for radius=[0 .65 .87]
        p=mu.poisson_density(theta,radius);
        verifyGreaterThanOrEqual(tc,min(p(:)),-1e-14);
        verifyLessThan(tc,abs(trapz(theta,p)-1),3e-7);
    end
    verifyLessThan(tc,max(abs(mu.poisson_density(theta,.8)- ...
        js_geronimus_poisson(theta+2*pi,a,.8))),2e-12);
end
end

function test_finite_cmv_cyclic_convergence_not_trace(tc)
% beta=+1 avoids exponentially close finite nodes near the gap atom.
a=.8;beta=1;r=.72;
theta=linspace(-pi,pi,301)';
target=js_geronimus_poisson(theta,a,r);
errs=zeros(3,1);
for j=1:3
    n=12*2^(j-1);
    U=full(geronimus_cmv(n,a,'cmv',beta));
    verifyLessThan(tc,norm(U'*U-eye(n),'fro'),2e-10);
    [z,w]=szego_quadrature(U);
    verifyLessThan(tc,abs(sum(w)-1),3e-12);
    verifyLessThan(tc,abs(sum(w.*z)-conj(a)),3e-10);
    errs(j)=max(abs(finite_poisson(theta,z,w,r)-target));
    trace_m1=trace(U)/n;
    verifyLessThan(tc,abs(trace_m1-mean(z)),3e-12);
    if n>=24
        verifyGreaterThan(tc,abs(trace_m1-conj(a)),.2);
    end
    S=geronimus_cmv(n,a,'cmv',beta,[],'spectral');
    verifyEqual(tc,S.finite.kind,'unknown'); % no explicit finite nodes
    verifyEqual(tc,S.underlying.kind,'analytic');
    verifyLessThan(tc,abs(S.underlying.atoms.weights-8/9),3e-15);
    verifyLessThan(tc,abs(S.underlying.poisson_density(.3,r)- ...
        js_geronimus_poisson(.3,a,r)),1e-14);
end
verifyLessThan(tc,errs(3),errs(2));
verifyLessThan(tc,errs(2),errs(1));
verifyLessThan(tc,errs(3),1e-4);
end

function test_rogers_szego_smooth_reference_control(tc)
q=.64;n=48;r=.7;
mu=szego_rogers_measure(q);
theta=linspace(-pi,pi,301)';
target=mu.circle_density(theta);
verifyGreaterThan(tc,min(target),0);
verifyLessThan(tc,abs(trapz(theta,target)-1),1e-8);
S=rogers_szego(n,q,'cmv',1,[],'spectral');
verifyEqual(tc,S.underlying.kind,'analytic');
U=full(rogers_szego(n,q));
[z,w]=szego_quadrature(U);
verifyLessThan(tc,abs(sum(w.*z)-sqrt(q)),4e-10);
verifyLessThan(tc,norm(U'*U-eye(n),'fro'),1e-10);
verifyLessThan(tc,abs(sum(w)-1),2e-12);
verifyGreaterThan(tc,min(finite_poisson(theta,z,w,r)),0);
end

function test_invalid_density_inputs(tc)
verifyError(tc,@()js_geronimus_density(0,1), ...
    'js_geronimus_density:Alpha');
verifyError(tc,@()js_geronimus_poisson(.1,1,.8), ...
    'js_geronimus_poisson:Alpha');
verifyError(tc,@()js_geronimus_poisson(.1,.2,1), ...
    'MATLAB:notLess');
end

function p=finite_poisson(theta,z,w,r)
% Equal to convolution with a discrete E1 cyclic measure, not trace.
d=theta(:)-angle(z(:)).';
K=(1-r^2)./(2*pi*(1-2*r*cos(d)+r^2));
p=K*w(:);
end
