function tests=test_tensor_fem_wave
%TEST_TENSOR_FEM_WAVE Exact 2D/3D tensor-product consistent-FE waves.
tests=functiontests(localfunctions);
end

function test_independent_2D_FEM_symbols_and_energy_metrics(tc)
for grid={[2 3],[3 2],[1 4],[3 3]}
    dims=grid{1};N=prod(dims);
    L=[1.4 .9];c=[1.25 .8];omega0=.63;
    h=.04*min(L./(dims+1)./c);
    [Mass,K]=reference_tensors(dims,L,c,omega0);
    physicalW2=sort(eig(full(K),full(Mass)));
    for method={'midpoint','gauss2','gauss3','verlet'}
        name=method{1};
        U=tensor_fem_wave(dims,L,c,omega0,h,name);
        P=tensor_fem_wave(dims,L,c,omega0,h,name,'pencil');
        S=tensor_fem_wave(dims,L,c,omega0,h,name,'spectral');
        verifyTrue(tc,isreal(U));
        verifyEqual(tc,size(U),[2*N 2*N]);
        verifyEqual(tc,size(P.M),[2*N 2*N]);
        verifyEqual(tc,S.dimension,2*N);
        verifyEqual(tc,S.spatial_dimension,2);
        verifyEqual(tc,S.grid_shape,dims);
        verifyEqual(tc,S.finite.kind,'exact');
        verifyLessThan(tc,norm(full(P.mass-Mass),'fro'),1e-12);
        verifyLessThan(tc,norm(full(P.stiffness-K),'fro'),1e-12);
        verifyLessThan(tc,norm(sort(S.frequencies.^2)-physicalW2), ...
            1e-10*max(1,max(physicalW2)));
        verifyLessThan(tc,norm(U'*U-eye(2*N),'fro'),2e-11*N);
        Jc=[zeros(N),eye(N);-eye(N),zeros(N)];
        verifyLessThan(tc,norm(U'*Jc*U-Jc,'fro'),2e-11*N);
        verifyLessThan(tc,norm(P.M'*Jc*P.M-Jc,'fro'), ...
            2e-9*N);
        verifyLessThan(tc,norm(P.M'*P.G*P.M-P.G,'fro')/ ...
            norm(P.G,'fro'),2e-9);
        verifyGreaterThan(tc,min(eig(P.G)),0);
        verifyLessThan(tc,norm(P.T'*P.T-P.G,'fro')/ ...
            norm(P.G,'fro'),2e-10);
        verifyLessThan(tc,norm(P.T*P.M/P.T-U,'fro'), ...
            2e-9*N);
        verifyLessThan(tc,norm(P.A'* (P.B\P.A)-P.B,'fro')/ ...
            norm(P.B,'fro'),2e-9);
        expectedGq=full(K);
        if strcmp(name,'verlet')
            expectedGq=full(K)-.25*h^2*full(K*(Mass\K));
        end
        expectedG=blkdiag(expectedGq,full(Mass\speye(N)));
        verifyLessThan(tc,norm(P.G-expectedG,'fro')/ ...
            norm(expectedG,'fro'),2e-9);
        eigenvalues=S.finite.eigenvalues();
        verifyEqual(tc,numel(eigenvalues),2*N);
        verifyLessThan(tc,max(arrayfun(@(x)min(abs(eig(U)-x)), ...
            eigenvalues)),3e-8);
        for p=-3:4
            verifyLessThan(tc,abs(trace(U^abs(p))/(2*N) ...
                -S.finite.moments(abs(p))),5e-10);
            verifyLessThan(tc,abs(mean(eigenvalues.^p) ...
                -S.finite.moments(p)),5e-10);
        end
    end
end
end

function test_exact_3D_FE_symbols_and_operator(tc)
for dims={[2 2 2],[1 2 3],[2 1 3]}
    grid=dims{1};N=prod(grid);
    L=[1.1 .8 1.4];c=[.9 1.1 .7];omega0=.2;
    [Mass,K]=reference_tensors(grid,L,c,omega0);
    h=.025;
    for method={'midpoint','gauss2','verlet'}
        name=method{1};
        U=tensor_fem_wave(grid,L,c,omega0,h,name);
        S=tensor_fem_wave(grid,L,c,omega0,h,name,'spectral');
        P=tensor_fem_wave(grid,L,c,omega0,h,name,'pencil');
        O=tensor_fem_wave(grid,L,c,omega0,h,name,'operator');
        verifyEqual(tc,S.spatial_dimension,3);
        verifyEqual(tc,S.dimension,2*N);
        verifyLessThan(tc,norm(sort(S.frequencies.^2) ...
            -sort(eig(full(K),full(Mass)))),2e-10);
        verifyLessThan(tc,norm(U'*U-eye(2*N),'fro'),5e-11*N);
        verifyLessThan(tc,norm(full(P.mass-Mass),'fro'),2e-12);
        verifyLessThan(tc,norm(full(P.stiffness-K),'fro'),2e-12);
        X=reshape(sin((1:4*N)*.31),2*N,2)+ ...
          1i*reshape(cos((1:4*N)*.17),2*N,2);
        verifyEqual(tc,O([],'size'),[2*N 2*N]);
        verifyLessThan(tc,norm(O(X,'notransp')-U*X,'fro'), ...
            3e-11*N);
        verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'), ...
            3e-11*N);
    end
end
end

function test_repeated_multidimensional_modes_and_phase_wrapping(tc)
n=4;grid=[n n];
S=tensor_fem_wave(grid,[1 1],[1 1],.2,.017, ...
        'midpoint','spectral');
F=reshape(S.frequencies,n,n);
verifyLessThan(tc,norm(F-F','fro'),1e-13);
verifyLessThan(tc,abs(F(2,3)-F(3,2)),1e-13);
R=tensor_fem_wave(grid,[1 1.3],[1 1],.2,.017, ...
        'midpoint','spectral');
A=reshape(R.frequencies,n,n);
verifyGreaterThan(tc,norm(A-A','fro'),.1);
for method={'midpoint','gauss2','gauss3'}
    S=tensor_fem_wave([3 4],[1 .8],[1 1.1],.1,2, ...
              method{1},'spectral');
    z=S.finite.eigenvalues();
    verifyLessThan(tc,max(abs(abs(z)-1)),2e-12);
    for k=[0 1 2 7 13]
        verifyLessThan(tc,abs(S.finite.moments(k) ...
            -mean(z.^k)),2e-12);
    end
end
end

function test_Verlet_mass_weighted_palindromic_QEP(tc)
for grid={[2 3],[2 2 2]}
    dims=grid{1};N=prod(dims);
    L=ones(1,numel(dims));c=L;
    S=tensor_fem_wave(dims,L,c,.3,.015,'verlet','spectral');
    Q=tensor_fem_wave(dims,L,c,.3,.015,'verlet','polynomial');
    verifyEqual(tc,Q.dimension,N);
    verifyLessThan(tc,norm(full(Q.A0-Q.A2),'fro'),1e-13);
    verifyLessThan(tc,norm(full(Q.A0-Q.A0'),'fro'),1e-13);
    verifyLessThan(tc,norm(full(Q.A1-Q.A1'),'fro'),1e-13);
    verifyGreaterThan(tc,min(eig(full(Q.A0))),0);
    W=[zeros(N),eye(N);-eye(N),-full(Q.A0\Q.A1)];
    z=eig(W);
    finite=S.finite.eigenvalues();
    verifyLessThan(tc,max(arrayfun(@(x)min(abs(z-x)),finite)),2e-7);
    P=tensor_fem_wave(dims,L,c,.3,.015,'verlet','pencil');
    Mass=full(P.mass);K=full(P.stiffness);
    eyeN=eye(N);
    actual=[eyeN-.5*.015^2*(Mass\K),.015*(Mass\eyeN); ...
        -.015*K+.25*.015^3*K*(Mass\K), ...
        eyeN-.5*.015^2*(K/Mass)];
    verifyLessThan(tc,norm(P.M-actual,'fro'),2e-8*N);
end
end

function test_Gauss_Pade_physical_step_and_order(tc)
% Compare the physical map against the exact Hamiltonian expm, not a
% reimplementation of modal wave phases.
dims=[1 2];N=prod(dims);
[Mass,K]=reference_tensors(dims,[1.1 .8],[1.3 .6],.5);
A=[zeros(N),full(Mass\speye(N));-full(K),zeros(N)];
for method={'midpoint','gauss2','gauss3'}
    name=method{1};
    h=.11;
    P=tensor_fem_wave(dims,[1.1 .8],[1.3 .6],.5,h,name,'pencil');
    coeff=[];
    switch name
        case 'midpoint',coeff=[1 .5];
        case 'gauss2',coeff=[1 .5 1/12];
        case 'gauss3',coeff=[1 .5 1/10 1/120];
    end
    Z=h*A;Pplus=zeros(2*N);Pminus=zeros(2*N);
    for j=1:numel(coeff)
        Pplus=Pplus+coeff(j)*Z^(j-1);
        Pminus=Pminus+coeff(j)*(-Z)^(j-1);
    end
    expected=Pminus\Pplus;
    verifyLessThan(tc,norm(P.M-expected,'fro'),5e-10);
    small=tensor_fem_wave(dims,[1.1 .8],[1.3 .6],.5,h/2, ...
        name,'pencil');
    err1=norm(P.M-expm(h*A),'fro');
    err2=norm(small.M-expm(.5*h*A),'fro');
    verifyGreaterThan(tc,err1,err2*4);
end
end

function test_large_fast_separable_operators(tc)
for grid={[48 37],[14 12 8],[40 35 7]}
    dims=grid{1};d=numel(dims);N=prod(dims);
    O=tensor_fem_wave(dims,ones(1,d),ones(1,d),.4,[], ...
        'gauss3','operator');
    verifyEqual(tc,O([],'size'),[2*N,2*N]);
    x=complex(sin(.7*(1:2*N)'),cos(.41*(1:2*N)'));
    y=O(x,'notransp');
    verifyLessThan(tc,norm(O(y,'ctransp')-x), ...
        4e-10*norm(x));
    verifyLessThan(tc,abs(norm(y)-norm(x)),4e-10*norm(x));
    verifyError(tc,@()O(ones(2*N-1,1),'notransp'), ...
        'um_tensor_fem_apply:Dimension');
    verifyError(tc,@()O(x,'oops'), ...
        'um_tensor_fem_apply:Flag');
end
end

function test_public_inventory_and_validation(tc)
[U,p]=tensor_fem_wave();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyTrue(tc,any(strcmp(p,'orthogonal')));
verifyTrue(tc,any(strcmp(p,'scalable')));
verifyEqual(tc,size(tensor_fem_wave([2 3])),[12 12]);
verifyError(tc,@()tensor_fem_wave(4),'tensor_fem_wave:Grid');
verifyError(tc,@()tensor_fem_wave([1 2 3 4]), ...
    'tensor_fem_wave:Grid');
verifyError(tc,@()tensor_fem_wave([2 0]), ...
    'tensor_fem_wave:Grid');
verifyError(tc,@()tensor_fem_wave([2 3],[1 2 3]), ...
    'tensor_fem_wave:Lengths');
verifyError(tc,@()tensor_fem_wave([2 3],[1 2],[0 1]), ...
    'tensor_fem_wave:Speeds');
verifyError(tc,@()tensor_fem_wave([2 3],[],[],[],[],'foo'), ...
    'tensor_fem_wave:Method');
verifyError(tc,@()tensor_fem_wave([2 3],[],[],[],[],'midpoint','polynomial'), ...
    'tensor_fem_wave:Polynomial');
verifyError(tc,@()tensor_fem_wave([12 12]), ...
    'tensor_fem_wave:ExplicitOrder');
S=tensor_fem_wave([2 3],[],[],[],[],'verlet','spectral');
verifyError(tc,@()tensor_fem_wave([2 3],[],[],[], ...
    2.1/max(S.frequencies),'verlet','operator'), ...
    'tensor_fem_wave:Stability');
end

function [Mass,K]=reference_tensors(grid,L,c,omega0)
% Independent direct 1D Galerkin matrices, then tensor products.
d=numel(grid);ones1=cell(1,d);mass=cell(1,d);stiff=cell(1,d);
for j=1:d
    n=grid(j);dx=L(j)/(n+1);
    D=zeros(n);G=zeros(n);
    for i=1:n
        D(i,i)=4*dx/6;
        G(i,i)=2/dx;
        if i<n
            D(i,i+1)=dx/6;D(i+1,i)=dx/6;
            G(i,i+1)=-1/dx;G(i+1,i)=-1/dx;
        end
    end
    mass{j}=D;stiff{j}=G;
    ones1{j}=eye(n);
end
Mass=1;
for j=d:-1:1,Mass=kron(Mass,mass{j});end
K=omega0^2*Mass;
for j=1:d
    B=1;
    for k=d:-1:1
        if k==j,F=stiff{k};else,F=mass{k};end
        B=kron(B,F);
    end
    K=K+c(j)^2*B;
end
Mass=sparse(Mass);K=sparse(K);
end
