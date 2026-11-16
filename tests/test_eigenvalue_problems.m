function tests=test_eigenvalue_problems
%TEST_EIGENVALUE_PROBLEMS Strict unitary pencils and nonlinear references.
tests=functiontests(localfunctions);
end

function test_pencil_invariants_and_property_contract(testCase)
cases={...
'mobius_unitary_pencil',{9,.7*exp(.3i)}; ...
'unitary_delay_pencil',{3,[2;3;4],[]}; ...
'finite_cmv',{7,[.2*ones(6,1);exp(.3i)]}; ...
'block_cmv_pencil',{4,2,[],[]}; ...
'floquet_factor_pencil',{5,.3,-.7}; ...
'cayley_hermitian_pencil',{9,.3,.7,2}; ...
'left_weighted_unitary_pencil',{9,1e7}};
for j=1:size(cases,1)
    name=cases{j,1};args=cases{j,2};
    [U,props,A,B]=feval(name,args{:});
    N=size(U,1);
    verifyEqual(testCase,size(U),[N,N],name);
    verifyEqual(testCase,size(A),[N,N],name);
    verifyEqual(testCase,size(B),[N,N],name);
    verifyTrue(testCase,any(strcmp('unitary',props)),name);
    verifyLessThan(testCase,norm(full(U'*U-speye(N)),'fro'), ...
        2e-9*max(1,N),name);
    verifyLessThan(testCase,norm(full(A*A'-B*B'),'fro')/ ...
        max(1,norm(full(B*B'),'fro')),1e-9,name);
    verifyLessThan(testCase,norm(full(B\A-U),'fro'), ...
        2e-9*max(1,N),name);
    verifyLessThan(testCase,max(abs(abs(eig(full(U)))-1)), ...
        5e-9,name);
    [empty,defaultProperties]=feval(name);
    verifyEmpty(testCase,empty,name);
    verifyEqual(testCase,props,defaultProperties,name);
end
end

function test_mobius_exact_and_poisson_limit(testCase)
for n=[1 2 7 18]
    for a=[0 .4 .75*exp(.4i)]
        [U,~,A,B]=mobius_unitary_pencil(n,a);
        omega=exp(2i*pi*(0:n-1)'/n);
        z=(omega-a)./(1-conj(a)*omega);
        verifyLessThan(testCase,abs(trace(U)-sum(z)),5e-11*n);
        verifyLessThan(testCase,max(abs(eig(full(U)).^n-1)) ...
            * double(a==0),1e-8);
        S=mobius_unitary_pencil(n,a,'spectral');
        verifyEqual(testCase,S.finite.kind,'exact');
        verifyEqual(testCase,S.dimension,n);
        verifyLessThan(testCase,abs(S.finite.moments(3)-mean(z.^3)),1e-12);
        verifyLessThan(testCase,abs(S.finite.moments(-2)-mean(z.^(-2))),1e-12);
        verifyLessThan(testCase,max(abs(S.limit.moments(0:5)-(-a).^(0:5))),1e-12);
        verifyLessThan(testCase,abs(integral(S.limit.density,-pi,pi)-1),2e-10);
        verifyLessThan(testCase,norm(full(A*A'-B*B'),'fro'),1e-10*n);
    end
end
end

function test_delay_secular_polynomial_and_scattering(testCase)
for lengths={[1;1],[2;3],[1;3;2]}
    ell=lengths{1};m=numel(ell);N=sum(ell);
    S=discrete_fourier(m);
    [U,~,A,B]=unitary_delay_pencil(m,ell,S);
    verifyEqual(testCase,size(U),[N N]);
    verifyTrue(testCase,issparse(U));
    z0=0; p0=det(z0*eye(N)-full(U));
    for z=[.3+.2i, -.7+.3i, 1.15+.1i]
        p=det(z*eye(N)-full(U))/p0;
        q=det(eye(m)-diag(z.^ell)*S);
        verifyLessThan(testCase,abs(p-q),2e-9*(1+abs(q)));
    end
    roots=eig(full(U));
    for k=1:N
        Z=diag(roots(k).^ell);
        verifyLessThan(testCase,abs(det(eye(m)-Z*S)),2e-8);
    end
    verifyLessThan(testCase,norm(full(A*A'-B*B'),'fro'),2e-12);
    R=unitary_delay_pencil(m,ell,S,'spectral');
    verifyEqual(testCase,R.dimension,N);
    verifyEqual(testCase,R.finite.kind,'unknown');
end
end

function test_scalar_cmv_factor_equivalence(testCase)
for n=[1 2 3 4 7 10]
    alpha=[.2*exp(1i*(1:n-1)');exp(.23i)];
    [U,~,A,B]=finite_cmv(n,alpha);
    C=finite_cmv(n,alpha);
    verifyLessThan(testCase,norm(full(U-C),'fro'),2e-11);
    verifyLessThan(testCase,norm(full(B\A-C),'fro'),2e-11);
    verifyTrue(testCase,issparse(U));
end
end

function test_block_cmv_matrix_coefficients(testCase)
q=2;
alpha=zeros(q,q,3);
alpha(:,:,1)=[.2,.1i;-.05,.1];
alpha(:,:,2)=[0,.27;-.11i,.08];
alpha(:,:,3)=[.18,.01;-.1,.2i];
beta=[0,1;-1,0];
for n=[1 2 3 4]
    [U,~,A,B]=block_cmv_pencil(n,q,alpha(:,:,1:max(0,n-1)),beta);
    verifyLessThan(testCase,norm(full(U'*U-eye(q*n)),'fro'),1e-11);
    verifyLessThan(testCase,norm(full(B\A-U),'fro'),1e-11);
    verifyTrue(testCase,issparse(U));
    spec=block_cmv_pencil(n,q, ...
        alpha(:,:,1:max(0,n-1)),beta,'spectral');
    verifyEqual(testCase,spec.dimension,n*q);
end
% The scalar specialization agrees with the canonical CMV implementation.
v=[.2+.1i;-.18;exp(.2i)];
C=block_cmv_pencil(3,1,reshape(v(1:2),1,1,2),v(3));
verifyLessThan(testCase,norm(full(C-finite_cmv(3,v)),'fro'),1e-11);
end

function test_cayley_and_floquet_compatibility(testCase)
for n=[1 2 5 9]
    [U,~,A,B]=cayley_hermitian_pencil(n,.1,.3,.8);
    V=cayley_jacobi(n,.1,.3,.8);
    verifyLessThan(testCase,norm(full(U-V),'fro'),3e-11*n);
    verifyLessThan(testCase,norm(full(A*A'-B*B'),'fro'),3e-11*n);
end
for n=[1 2 3 6]
    [U,~,A,B]=floquet_factor_pencil(n,.43,-.61);
    V=split_step_walk(n,.43,-.61);
    verifyLessThan(testCase,norm(full(U-V),'fro'),5e-12*n);
    verifyLessThan(testCase,norm(full(B\A-V),'fro'),5e-12*n);
end
end

function test_weighted_ill_scaled_pencil(testCase)
for n=[1 2 9]
    for kappa=[1,1e4,1e8]
        [U,~,A,B]=left_weighted_unitary_pencil(n,kappa);
        verifyEqual(testCase,U,cycle_shift(n));
        verifyLessThan(testCase,norm(full(A*A'-B*B'),'fro')/ ...
            max(1,norm(full(B*B'),'fro')),2e-12);
        verifyLessThan(testCase,norm(full(B\A-U),'fro'),2e-12);
    end
end
end

function test_invalid_pencil_parameters(testCase)
verifyError(testCase,@()mobius_unitary_pencil(5,1), ...
    'mobius_unitary_pencil:Parameter');
verifyError(testCase,@()unitary_delay_pencil(2,[2;0]), ...
    'MATLAB:expectedPositive');
verifyError(testCase,@()unitary_delay_pencil(2,[2;3],2*eye(2)), ...
    'unitary_delay_pencil:Scattering');
verifyError(testCase,@()unitary_delay_pencil(2,[1;2;3]), ...
    'unitary_delay_pencil:Delays');
verifyError(testCase,@()finite_cmv(3,[.2;1;1]), ...
    'js_core_product:Alpha');
verifyError(testCase,@()block_cmv_pencil(3,2,zeros(2,2,2),2*eye(2)), ...
    'block_cmv_pencil:Beta');
a=zeros(2,2,2);a(:,:,1)=eye(2);
verifyError(testCase,@()block_cmv_pencil(3,2,a), ...
    'block_cmv_pencil:Alpha');
verifyError(testCase,@()cayley_hermitian_pencil(4,0,0,0), ...
    'MATLAB:expectedPositive');
verifyError(testCase,@()left_weighted_unitary_pencil(4,.5), ...
    'MATLAB:expectedGreaterThanOrEqual');
end

function test_eigenproblem_registry_integrity(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
catalog=jsondecode(fileread(fullfile(root,'metadata','catalog.json')));
reg=jsondecode(fileread(fullfile(root,'metadata','eigenvalue_problems.json')));
verifyEqual(testCase,numel(reg.families),14);
classes={'strict_unitary_pencil','unitary_on_contour', ...
         'circle_spectrum_only','reciprocal_only'};
ids={catalog.matrix_families.id};
implemented=0;
for j=1:numel(reg.families)
    e=reg.families(j);
    verifyTrue(testCase,any(strcmp(e.classification,classes)));
    if strcmp(e.status,'implemented') || strcmp(e.status,'existing_generator')
        implemented=implemented+1;
        verifyTrue(testCase,any(strcmp(e.id,ids)));
        verifyTrue(testCase,isfile(fullfile(root,[e.id '.m'])));
        category=catalog.matrix_families(strcmp(ids,e.id)).category;
        if strcmp(e.status,'existing_generator')
            verifyEqual(testCase,category,'opuc');
        else
            verifyEqual(testCase,category,'eigenvalue_problems');
        end
    else
        verifyFalse(testCase,any(strcmp(e.id,ids)), ...
            'External NEPs must not be labeled as unitary matrix generators.');
    end
end
verifyEqual(testCase,implemented,7);
end
