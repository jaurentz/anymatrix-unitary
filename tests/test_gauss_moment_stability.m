function tests=test_gauss_moment_stability
%TEST_GAUSS_MOMENT_STABILITY Independent normalized moments of GH/GL.
tests=functiontests(localfunctions);
end

function test_hermite_analytic_moments_and_exactness_order(tc)
exact=[1;0;.5;0;.75;0;15/8;0;105/16];
for n=[1 2 3 8 32 96 192 384]
    S=gauss_hermite_transform(n,'quadrature');
    M=S.moment_diagnostics;
    degree=min(8,2*n-1);
    verifyEqual(tc,M.moment_measure,'hermite');
    verifyEqual(tc,M.moment_degree_exact_in_arithmetic,2*n-1);
    verifyEqual(tc,M.moment_orders,(0:degree)');
    verifyLessThan(tc,norm(M.exact_moments-exact(1:degree+1)), ...
        3e-14);
    verifyTrue(tc,M.moment_finite);
    verifyLessThan(tc,M.maximum_scaled_moment_error, ...
        2e-8);
    verifyLessThan(tc,abs(M.quadrature_moments(1)-1), ...
        3e-9);
    verifyLessThan(tc,norm(M.quadrature_moments ...
        -M.exact_moments),1e-7);
end
end

function test_gamma_analytic_moments_extreme_alpha(tc)
for alpha=[-.999999 -.999 -.7 0 .5 3 16]
    for n=[1 2 3 16 64 192]
        S=gauss_laguerre_transform(n,alpha,'quadrature');
        M=S.moment_diagnostics;
        degree=min(8,2*n-1);
        exact=ones(degree+1,1);
        for k=1:degree
            exact(k+1)=prod(alpha+(1:k));
        end
        verifyEqual(tc,M.moment_orders,(0:degree)');
        verifyEqual(tc,M.moment_degree_exact_in_arithmetic,2*n-1);
        verifyEqual(tc,M.moment_measure,'laguerre');
        verifyLessThan(tc,norm(M.exact_moments-exact), ...
            2e-12*max(1,norm(exact)));
        verifyTrue(tc,M.moment_finite);
        verifyLessThan(tc,M.maximum_scaled_moment_error, ...
            2e-6);
        verifyLessThan(tc,abs(M.quadrature_moments(1)-1), ...
            3e-7);
    end
end
end

function test_moments_are_quadrature_not_unitary_spectrum(tc)
n=10;
H=gauss_hermite_transform(n,'quadrature');
L=gauss_laguerre_transform(n,.37,'quadrature');
S1=gauss_hermite_transform(n,'spectral');
S2=gauss_laguerre_transform(n,.37,'spectral');
verifyEqual(tc,S1.finite.kind,'unknown');
verifyEqual(tc,S2.finite.kind,'unknown');
verifyTrue(tc,all(isreal(H.nodes)));
verifyTrue(tc,all(isreal(L.nodes)));
verifyGreaterThan(tc,min(L.nodes),0);
verifyLessThan(tc,H.moment_diagnostics.maximum_scaled_moment_error,1e-9);
verifyLessThan(tc,L.moment_diagnostics.maximum_scaled_moment_error,1e-9);
verifyError(tc,@()um_gauss_moment_diagnostics(1,0,'other'), ...
    'um_gauss_moment_diagnostics:Family');
verifyError(tc,@()um_gauss_moment_diagnostics(1,0,'laguerre',-1), ...
    'um_gauss_moment_diagnostics:Alpha');
end
