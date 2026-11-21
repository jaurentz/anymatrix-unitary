function Y=um_mesh_spin_dirac_cayley(X,flag,N,steps,plus,minus,plusMatrix,minusMatrix)
%UM_MESH_SPIN_DIRAC_CAYLEY Sparse rational mesh spin Dirac evolution.
%   For Hermitian H and h=t/steps, C=(I+i*h*H/2)\(I-i*h*H/2)
%   is unitary in mathematical arithmetic and C^steps approximates
%   exp(-i*t*H) to second order at fixed H. Because H is REAL
%   symmetric here, C is complex symmetric: C.'=C. The ordinary
%   transpose action therefore equals forward. The adjoint is C^-1.
%   Floating point sparse LU solves and finite-step rational error
%   must not be described as an exact exponential or exact isometry.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_mesh_spin_dirac_cayley:Flag','Flag must be scalar text.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[N N];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_mesh_spin_dirac_cayley:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_mesh_spin_dirac_cayley:Dimension', ...
        'Input must have 4*nF rows.');
end
Y=full(X);
for k=1:steps
    if strcmp(flag,'ctransp')
        Y=minus\(plusMatrix*Y);
    else
        Y=plus\(minusMatrix*Y);
    end
end
end
