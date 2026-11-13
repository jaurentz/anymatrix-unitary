function tests=test_spectral_reference
%TEST_SPECTRAL_REFERENCE Analytical reference laws and API regression.
tests=functiontests(localfunctions);
end

function test_all_public_generators_support_mode(testCase)
cases={ ...
 'cycle_shift',{7}; 'weighted_shift',{7,[]}; ...
 'diagonal_phase',{7,[]}; 'phase_circulant',{7,[]}; ...
 'cycle_propagator',{7,[]}; 'orthogonal_blocks',{7,[],[]}; ...
 'finite_cmv',{7,[]}; 'opuc_matrix',{7,[],[],[],[]}; ...
 'jacobi_szego',{7,[],[],[],[],[]}; ...
 'bernstein_szego',{7,[],[],[],[]}; ...
 'circular_beta',{7,[],[]}; ...
 'circular_jacobi',{7,[],[],[],[],[]}; ...
 'rogers_szego',{7,[],[],[],[]}; ...
 'geronimus_cmv',{7,[],[],[],[]}; ...
 'periodic_cmv',{7,[],[],[],[]}; ...
 'bernstein_szego_finite',{7,[],[],[],[]}; ...
 'haar_unitary',{7,[]}; 'haar_orthogonal',{7,[]}; ...
 'haar_wavelet',{8,[]}; 'orthogonal_dct',{8,[]}; ...
 'discrete_hartley',{8}; 'cazac_circulant',{7,[]}; ...
 'fresnel_transform',{8}; 'affine_fourier',{8,[],[]}; ...
 'finite_zak',{8,[]}; 'time_frequency_shift',{8,[],[]}; ...
 'periodic_hilbert',{8,[],[]}; 'discrete_fourier',{8}; ...
 'discrete_sine',{8,[]}; 'walsh_hadamard',{8,[]}; ...
 'dpss_transform',{8,[]}; 'graph_fourier',{8,[]}; ...
 'split_step_walk',{7,[],[]}; 'dirac_automaton',{7,[]}; ...
 'grover_operator',{7,[]}; 'quantum_baker',{8,[]}; ...
 'kicked_harper',{8,[],[]}; 'kicked_top',{8,[],[]}; ...
 'floquet_ising',{4,[],[],[]}; ...
 'selfdual_kicked_ising',{4,[]}; ...
 'driven_kitaev',{7,[]}; ...
 'kicked_rotor',{8,[],[],[]}; ...
 'dual_unitary_circuit',{4,[],[],[],[]}; ...
 'xxz_floquet',{4,[],[],[]}; ...
 'kicked_bose_hubbard',{7,[],[],[],[],[]}; ...
 'jaynes_cummings',{7,[],[],[],[]}; ...
 'chalker_coddington',{2,3,[],[]}; ...
 'cayley_jacobi',{7,[],[],[]}; ...
 'polar_unitary',{7,[]}; 'julia_unitary',{7,[]}};
verifyEqual(testCase,size(cases,1),50);
for j=1:size(cases,1)
    name=cases{j,1}; inputs=cases{j,2};
    S=feval(name,inputs{:},'spectral');
    verifyTrue(testCase,isstruct(S),name);
    verifyEqual(testCase,S.family,name);
    verifyGreaterThan(testCase,S.dimension,0,name);
    verifyTrue(testCase,isfield(S,'finite')&&isfield(S,'limit') ...
        &&isfield(S,'ensemble')&&isfield(S,'underlying'),name);
end
end

function test_finite_spectra_match_independent_eig(testCase)
cases={ ...
 'cycle_shift',{17}; ...
 'weighted_shift',{11,.2*(0:10)'}; ...
 'diagonal_phase',{13,linspace(-2,2,13)'}; ...
 'phase_circulant',{13,linspace(-2,2,13)'}; ...
 'orthogonal_blocks',{9,[.1;.5;.9;1.5],-1}; ...
 'cycle_propagator',{11,.8}; ...
 'dirac_automaton',{11,.45}; ...
 'grover_operator',{17,3}; ...
 'discrete_fourier',{11}; ...
 'discrete_fourier',{12}; ...
 'walsh_hadamard',{16,'natural'}; ...
 'time_frequency_shift',{12,3,4}; ...
 'time_frequency_shift',{11,5,6}; ...
 'jaynes_cummings',{6,.8,.35,.7,1.2}};
for j=1:size(cases,1)
    name=cases{j,1}; inputs=cases{j,2};
    A=feval(name,inputs{:});
    S=feval(name,inputs{:},'spectral');
    verifyEqual(testCase,S.dimension,size(A,1),name);
    verifyEqual(testCase,S.finite.kind,'exact',name);
    z=S.finite.eigenvalues();
    actual=eig(full(A));
    verifyEqual(testCase,numel(z),numel(actual),name);
    for k=-3:5
        verifyLessThan(testCase,abs(mean(actual.^k)-S.finite.moments(k)), ...
            3e-9,name);
    end
    verifyLessThan(testCase,max(arrayfun(@(a)min(abs(a-z)),actual)), ...
        3e-9,name);
    verifyLessThan(testCase,max(arrayfun(@(a)min(abs(a-actual)),z)), ...
        3e-9,name);
end
end

function test_dft_all_four_residue_classes(testCase)
for n=1:21
    S=discrete_fourier(n,'spectral');
    z=S.finite.eigenvalues();
    F=discrete_fourier(n);
    verifyEqual(testCase,numel(z),n);
    verifyLessThan(testCase,abs(sum(z)-trace(F)),1e-10);
    verifyLessThan(testCase,abs(sum(z.^2)-trace(F^2)),1e-10);
    verifyLessThan(testCase,abs(sum(S.finite.atoms.weights)-1),1e-14);
end
end

function test_dirac_band_density_and_limits(testCase)
m=.45; S=dirac_automaton(30,m,'spectral');
a=asin(m);
verifyEqual(testCase,S.dimension,60);
verifyEqual(testCase,S.limit.kind,'analytic');
verifyLessThan(testCase,norm(S.limit.support- ...
    [-pi+a,-a;a,pi-a],'fro'),1e-12);
verifyEqual(testCase,S.limit.density(0),0);
verifyEqual(testCase,S.limit.density(pi),0);
verifyGreaterThan(testCase,S.limit.density(pi/2),0);
mass=0;
for j=1:2
    support=S.limit.support(j,:);
    mass=mass+integral(S.limit.density,support(1),support(2), ...
        'RelTol',1e-6,'AbsTol',1e-8);
end
verifyLessThan(testCase,abs(mass-1),3e-5);
S0=dirac_automaton(8,0,'spectral');
verifyLessThan(testCase,abs(S0.limit.density(.3)-1/(2*pi)),1e-14);
S1=dirac_automaton(8,1,'spectral');
verifyEqual(testCase,S1.limit.kind,'atomic');
verifyEqual(testCase,sum(S1.limit.atoms.weights),1);
end

function test_wrapped_arcsine_and_moments(testCase)
for t=[.35 2.6]
    S=cycle_propagator(12,t,'spectral');
    theta=[-2.8 -.75 0 .6 2.1];
    verifyLessThan(testCase,norm(S.limit.density(theta) ...
        -S.limit.density(theta+2*pi)),1e-10);
    for k=0:4
        verifyLessThan(testCase,abs(S.limit.moments(k) ...
            -besselj(0,2*t*k)),1e-12);
        samples=exp(-2i*t*k*cos(2*pi*(0:4095)'/4096));
        verifyLessThan(testCase,abs(mean(samples)-S.limit.moments(k)),1e-11);
    end
end
S=cycle_propagator(12,0,'spectral');
verifyEqual(testCase,S.limit.kind,'atomic');
end

function test_grover_multiplicities(testCase)
n=31;m=7;
S=grover_operator(n,m,'spectral');
verifyEqual(testCase,S.finite.atoms.weights(:), ...
    [1;1;m-1;n-m-1]/n,'AbsTol',1e-14);
z=S.finite.eigenvalues();
verifyEqual(testCase,sum(abs(z-1)<1e-12),m-1);
verifyEqual(testCase,sum(abs(z+1)<1e-12),n-m-1);
end

function test_underlying_opuc_is_not_trace(testCase)
S=jacobi_szego(16,0,0,'cmv',1,[],'spectral');
verifyEqual(testCase,S.finite.kind,'unknown');
verifyEqual(testCase,S.underlying.kind,'analytic');
mu=szego_jacobi_measure(0,0);
theta=linspace(-pi,pi,201);
verifyLessThan(testCase,max(abs(S.underlying.density(theta) ...
    -mu.circle_density(theta))),1e-13);
verifyLessThan(testCase,abs(integral(S.underlying.density,-pi,pi)-1),1e-10);
B=bernstein_szego(10,.2i,'cmv',1,[],'spectral');
R=rogers_szego(10,.5,'cmv',1,[],'spectral');
J=circular_jacobi(10,.3,.2,'cmv',1,[],'spectral');
verifyEqual(testCase,B.underlying.kind,'analytic');
verifyEqual(testCase,R.underlying.kind,'analytic');
verifyEqual(testCase,J.underlying.kind,'analytic');
end

function test_ensemble_rng_and_lazy_eigenvalues(testCase)
state=rng;
S=haar_unitary(100000,[],'spectral');
verifyEqual(testCase,rng,state);
verifyEqual(testCase,S.ensemble.kind,'exact_expectation');
verifyEqual(testCase,S.ensemble.density([0 pi]),ones(1,2)/(2*pi));
T=circular_beta(100000,2,123,'spectral');
verifyEqual(testCase,rng,state);
verifyEqual(testCase,T.finite.kind,'unknown');
Sbig=dirac_automaton(1e8,.4,'spectral');
verifyEqual(testCase,Sbig.dimension,2e8);
verifyTrue(testCase,isa(Sbig.finite.eigenvalues,'function_handle'));
end

function test_legacy_properties_and_operator(testCase)
[U,p]=dirac_automaton(6,.4);
verifyEqual(testCase,size(U),[12 12]);
verifyTrue(testCase,ismember('unitary',p));
[~,p0]=dirac_automaton();
verifyEqual(testCase,p,p0);
verifyTrue(testCase,isa(dirac_automaton(6,.4,'operator'), ...
    'function_handle'));
S=dirac_automaton(6,.4,'spectral');
verifyEqual(testCase,S.dimension,12);
[Q,q]=orthogonal_blocks(5,[.2;.3],-1);
verifyEqual(testCase,size(Q),[5 5]);
verifyTrue(testCase,ismember('orthogonal',q));
end
