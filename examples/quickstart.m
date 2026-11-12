%QUICKSTART Minimal Jacobi--Szego example (run from repository root).
mu = szego_jacobi_measure(0,0);
alpha = szego_verblunsky(mu,15);
U = szego_matrix(mu,16,'cmv',1);
H = szego_matrix(mu,16,'hessenberg',1);
[nodes,weights] = szego_quadrature(U);
disp(['Unitarity residual: ',num2str(norm(full(U'*U-eye(16)),'fro'))]);
disp(['Coefficient count: ',num2str(numel(alpha))]);
disp(['Quadrature weight sum: ',num2str(sum(weights))]);
assert(norm(full(U'*U-eye(16)),'fro') < 1e-12)
assert(abs(sum(weights)-1) < 1e-12)
assert(abs(trace(U)-trace(H)) < 1e-12)
