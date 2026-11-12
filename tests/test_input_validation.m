function tests=test_input_validation
%TEST_INPUT_VALIDATION Boundary values and errors.
tests=functiontests(localfunctions);
end

function test_wrong_dimensions(testCase)
verifyError(testCase,@()finite_cmv(5,[0;0;1]),'finite_cmv:DimensionMismatch');
verifyError(testCase,@()opuc_matrix(5,[0;0],1),'opuc_matrix:AlphaSize');
verifyError(testCase,@()diagonal_phase(5,[1;2]),'diagonal_phase:DimensionMismatch');
end

function test_unimodular_boundary(testCase)
mu=szego_jacobi_measure(.3,.4);
verifyError(testCase,@()szego_matrix(mu,5,'cmv',0.8),'js_core_product:Beta');
verifyError(testCase,@()finite_cmv(3,[0;0;0.8]),'js_core_product:Beta');
end

function test_default_construction(testCase)
for n=[1 2 3 7]
    U=jacobi_szego(n);
    verifyTrue(testCase,isreal(U));
    verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),2e-12*n);
end
end
