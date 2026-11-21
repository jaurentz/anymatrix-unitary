function tests=test_gauss_sturm_nodes
%TEST_GAUSS_STURM_NODES Eigenvector-free Gaussian nodes and weights (#40).
tests=functiontests(localfunctions);
end

function test_hermite_against_full_dense_gauss_basis(tc)
for n=[1 2 3 8 32 96 192]
    Z=gauss_hermite_transform(n,'nodes_weights');
    F=gauss_hermite_transform(n,'quadrature');
    verifyTrue(tc,Z.quadrature_only);
    verifyFalse(tc,Z.dense_transform_constructed);
    verifyTrue(tc,Z.sturm_diagnostics.converged);
    verifyFalse(tc,Z.sturm_diagnostics.dense_eigenvectors_constructed);
    verifyEqual(tc,numel(Z.nodes),n);
    verifyLessThan(tc,max(abs(Z.nodes-F.nodes)), ...
        5e-11*max(1,max(abs(F.nodes))));
    verifyLessThan(tc,abs(sum(Z.normalized_weights)-1),3e-8);
    verifyTrue(tc,all(isfinite(Z.log_normalized_weights)));
    verifyLessThan(tc,Z.moment_diagnostics.maximum_scaled_moment_error, ...
        5e-8);
    reliable=F.normalized_weights>1e-9;
    verifyLessThan(tc,max(abs(Z.log_normalized_weights(reliable) ...
        -F.log_normalized_weights(reliable))),5e-7);
    verifyLessThan(tc,norm(Z.log_function_scaling(reliable) ...
        -F.log_function_scaling(reliable)),5e-7);
    if n>1
        verifyTrue(tc,all(diff(Z.nodes)>0));
        verifyTrue(tc,all(Z.sturm_diagnostics.bracket_widths>=0));
    end
end
end

function test_generalized_laguerre_near_minus_one_and_tail(tc)
for alpha=[-.999999 -.7 0 .5 3 16]
    for n=[1 2 5 24 96 192]
        Z=gauss_laguerre_transform(n,alpha,'nodes_weights');
        D=gauss_laguerre_transform(n,alpha,'quadrature');
        verifyEqual(tc,Z.family,'gauss_laguerre_transform');
        verifyEqual(tc,Z.alpha,alpha);
        verifyTrue(tc,Z.quadrature_only);
        verifyTrue(tc,Z.sturm_diagnostics.converged);
        verifyGreaterThan(tc,min(Z.nodes),0);
        verifyLessThan(tc,max(abs(Z.nodes-D.nodes)), ...
            5e-11*max(1,max(abs(D.nodes))));
        verifyLessThan(tc,abs(sum(Z.normalized_weights)-1), ...
            6e-7);
        verifyTrue(tc,all(isfinite(Z.log_normalized_weights)));
        verifyTrue(tc,Z.moment_diagnostics.moment_finite);
        verifyLessThan(tc,Z.moment_diagnostics.maximum_scaled_moment_error, ...
            3e-6);
        reliable=D.normalized_weights>1e-8;
        verifyLessThan(tc,max(abs(Z.log_normalized_weights(reliable) ...
            -D.log_normalized_weights(reliable))),1e-4);
    end
end
end

function test_large_order_sturm_no_dense_matrix(tc)
% The numerical algorithm stores only O(n) vectors, no O(n^2) Q.
n=384;
for k=1:2
    if k==1
        Z=gauss_hermite_transform(n,'nodes_weights');
    else
        Z=gauss_laguerre_transform(n,.37,'nodes_weights');
    end
    verifyTrue(tc,Z.quadrature_only);
    verifyEqual(tc,numel(Z.sturm_diagnostics.iterations),n);
    verifyLessThan(tc,Z.sturm_diagnostics.maximum_bracket_width, ...
        1e-10*max(1,max(abs(Z.nodes))));
    verifyLessThan(tc,Z.moment_diagnostics.maximum_scaled_moment_error, ...
        5e-6);
    verifyLessThan(tc,abs(sum(Z.normalized_weights)-1),8e-7);
    verifyLessThan(tc,min(Z.log_normalized_weights),-700);
end
end

function test_arbitrary_irreducible_jacobi_and_invalid_inputs(tc)
n=51;
d=linspace(-2,3,n)';
e=.1+abs(sin((1:n-1)'))/2;
[x,status]=um_gauss_sturm_nodes(d,e);
J=diag(d)+diag(e,1)+diag(e,-1);
reference=sort(eig(J));
verifyLessThan(tc,max(abs(x-reference)),2e-11);
verifyTrue(tc,status.converged);
verifyLessThan(tc,max(status.bracket_widths),2e-12);
verifyError(tc,@()um_gauss_sturm_nodes([0;1],[0]), ...
    'um_gauss_sturm_nodes:Jacobi');
verifyError(tc,@()um_gauss_sturm_nodes([0;1],[NaN]), ...
    'um_gauss_sturm_nodes:Jacobi');
verifyError(tc,@()gauss_hermite_transform(4,'operator'), ...
    'gauss_hermite_transform:Mode');
verifyError(tc,@()gauss_laguerre_transform(4,.5,'operator'), ...
    'gauss_laguerre_transform:Mode');
end
