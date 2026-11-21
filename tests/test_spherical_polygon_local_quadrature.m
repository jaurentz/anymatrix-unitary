function tests=test_spherical_polygon_local_quadrature
%TEST_SPHERICAL_POLYGON_LOCAL_QUADRATURE Local versus global cubature (#144).
tests=functiontests(localfunctions);
end

function test_local_refinement_and_global_independent_q(tc)
shapes={[-.5 -.4;.57 -.38;.58 .45;-.53 .48], ...
        [-.71 -.49;.65 -.48;.12 -.06;.69 .53;-.6 .45]};
for k=1:numel(shapes)
    V=unit_north(shapes{k});
    for L=[1 2 4 6]
        n=L^2;
        A=spherical_polygon_slepian(L,V,1e-9,'adaptive_parameters');
        B=spherical_polygon_slepian(L,V,1e-9,'parameters');
        E=A.local_quadrature_diagnostics;
        verifyTrue(tc,A.quadrature_converged);
        verifyFalse(tc,A.quadrature_certified);
        verifyTrue(tc,E.converged);
        verifyFalse(tc,E.certified);
        verifyTrue(tc,all(E.triangle_converged));
        verifyEqual(tc,numel(E.triangle_orders),size(A.triangles,1));
        verifyTrue(tc,all(E.triangle_orders>=16));
        verifyTrue(tc,all(E.triangle_orders<=64));
        verifyTrue(tc,all(mod(E.triangle_orders,2)==0));
        verifyTrue(tc,all(E.triangle_levels>=2));
        verifyLessThan(tc,E.sum_exact_triangle_areas_error,3e-12);
        verifyLessThan(tc,abs(sum(E.triangle_exact_areas)-A.geometric_area), ...
            3e-12);
        verifyLessThan(tc,E.total_matrix_difference_upper_sum, ...
            1e-9*E.global_scale);
        verifyLessThan(tc,E.total_area_error,1e-9*E.global_scale);
        verifyLessThan(tc,E.global_addition_trace_error, ...
            1e-9*E.global_scale);
        verifyLessThan(tc,norm(A.concentration-B.concentration,'fro'), ...
            2e-8*n);
        verifyLessThan(tc,norm(A.basis'*A.basis-eye(n),'fro'), ...
            1e-10*n);
        verifyLessThan(tc,norm(A.basis'*A.concentration*A.basis ...
            -diag(A.concentration_eigenvalues),'fro'),4e-9*n);
        verifyLessThan(tc,norm(A.concentration ...
            -spherical_polygon_slepian(L,V,1e-9, ...
            'adaptive_concentration'),'fro'),4e-13*n);
    end
end
end

function test_nonconverged_status_and_strict_failure(tc)
V=unit_north([-.47 -.39;.53 -.35;.43 .56;-.54 .5]);
D=spherical_polygon_slepian(4,V,1e-24,'adaptive_diagnostics');
verifyFalse(tc,D.converged);
verifyFalse(tc,D.certified);
verifyEqual(tc,D.bandlimit,4);
verifyEqual(tc,size(D.concentration),[16 16]);
verifyGreaterThan(tc,D.total_matrix_difference_upper_sum,1e-24);
verifyTrue(tc,all(D.triangle_orders<=64));
verifyTrue(tc,all(isfinite(D.triangle_area_errors)));
verifyError(tc,@()spherical_polygon_slepian(4,V,1e-24, ...
    'adaptive_parameters'),'spherical_polygon_slepian:Quadrature');
verifyError(tc,@()spherical_polygon_slepian(4,V,1e-24, ...
    'adaptive_concentration'),'spherical_polygon_slepian:Quadrature');
end

function test_rotation_geometry_and_local_addition_theorem(tc)
V=unit_north([-.5 -.4;.62 -.41;.27 .06;.65 .48;-.52 .56]);
G=um_spherical_polygon_geometry(V);
L=5;
A=spherical_polygon_slepian(L,V,1e-9,'adaptive_parameters');
verifyLessThan(tc,norm(A.concentration-A.concentration','fro'), ...
    3e-12*L^2);
verifyLessThan(tc,abs(trace(A.concentration) ...
    -L^2*G.area/(4*pi)),3e-9*L^2);
verifyLessThan(tc,abs(A.numerical_area-G.area),3e-9);
verifyLessThan(tc,abs(A.local_quadrature_diagnostics ...
    .total_matrix_difference_upper_sum ...
    -sum(A.local_quadrature_diagnostics ...
    .triangle_frobenius_differences)),2e-14);
end

function V=unit_north(P)
V=[P ones(size(P,1),1)];
V=V./sqrt(sum(V.^2,2));
end
