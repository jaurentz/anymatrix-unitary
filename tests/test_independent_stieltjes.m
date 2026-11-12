function tests = test_independent_stieltjes
%TEST_INDEPENDENT_STIELTJES Jacobi-to-Szego regression entirely in MATLAB.
%   Gauss--Jacobi nodes and weights are constructed by Golub--Welsch.
%   A distinct discrete Stieltjes orthogonalization reconstructs their
%   recurrence; inverse Geronimus is compared against explicit OPUC data.
tests = functiontests(localfunctions);
end

function test_numerical_quadrature_to_verblunsky(testCase)
parameters = [0 0; -0.5 -0.5; 0.5 0.5; 0.3 1.6; -0.6 0.4; 2.4 -0.7];
for j = 1:size(parameters,1)
    a = parameters(j,1);
    b = parameters(j,2);
    [diagReference,subReference] = js_jacobi_recurrence(a,b,80);
    J = diag(diagReference) + diag(subReference,1) + diag(subReference,-1);
    [V,D] = eig(J);
    nodes = real(diag(D));
    weights = abs(V(1,:).').^2;
    weights = weights / sum(weights);
    verifyLessThan(testCase,abs(sum(weights)-1),1e-13);
    verifyGreaterThan(testCase,min(weights),0);

    [d,s] = local_stieltjes(nodes,weights,13);
    alpha = js_inverse_geronimus(d,s,24);
    expected = szego_verblunsky(szego_jacobi_measure(a,b),24);
    verifyLessThan(testCase,norm(alpha-expected,inf),1e-8);

    % The quadrature measure independently verifies circle Fourier moments.
    reference = js_jacobi_cosine_moments(a,b,6);
    t = max(-1,min(1,nodes/2));
    for k = 0:6
        moment = sum(weights.*cos(k*acos(t)));
        verifyLessThan(testCase,abs(moment-reference(k+1)),1e-8);
    end
end
end

function [d,s] = local_stieltjes(nodes,weights,m)
% Recover recurrence coefficients from a positive *discrete* measure.
nodes = nodes(:); weights = weights(:);
d = zeros(m,1);
s = zeros(m-1,1);
prev = zeros(size(nodes));
curr = ones(size(nodes));
aPrev = 0;
for k = 1:m
    d(k) = sum(weights.*nodes.*curr.^2);
    next = (nodes-d(k)).*curr - aPrev*prev;
    aNext = sqrt(sum(weights.*next.^2));
    if k<m, s(k)=aNext; end
    prev = curr;
    curr = next/aNext;
    aPrev = aNext;
end
end
