function tests=test_stiefel_parallel_transport
%TEST_STIEFEL_PARALLEL_TRANSPORT Canonical tangent parallel transport (#55).
tests=functiontests(localfunctions);
end

function test_orthogonality_tangent_frame_and_velocity(tc)
cases=[3 1;4 2;5 3;4 4;7 2];
for h=1:size(cases,1)
    n=cases(h,1);p=cases(h,2);
    A=zeros(p);
    if p>1
        A(1,2)=-.37;A(2,1)=.37;
    end
    B=zeros(n-p,p);
    if n>p
        B(1,1)=.45;
        if p>1,B(1,2)=.19;end
        if n-p>1,B(2,p)=-.23;end
    end
    t=.71;
    G=stiefel_parallel_transport(n,p,A,B,t,'model');
    d=n*p-p*(p+1)/2;
    verifyEqual(tc,size(G.U),[d d]);
    verifyLessThan(tc,norm(G.U'*G.U-eye(d),'fro'), ... 
        2e-11*max(1,d));
    verifyLessThan(tc,norm(G.connection_generator ...
        +G.connection_generator','fro'),2e-11);
    verifyLessThan(tc,norm(G.Yt'*G.Yt-eye(p),'fro'),2e-11);
    % Build tangent columns and check their canonical metric Gram matrices
    % at the start and at the moving endpoint independently.
    Z0=zeros(n*p,d);
    Zt=zeros(n*p,d);
    for k=1:d
        F=horizontal_core(n,G.basis_row(k),G.basis_col(k));
        Z0(:,k)=reshape(F(:,1:p),[],1);
        Fend=G.moving_frame*F(:,1:p);
        Zt(:,k)=reshape(Fend,[],1);
        verifyLessThan(tc,norm(G.Yt'*Fend+Fend'*G.Yt,'fro'), ...
            2e-11);
    end
    metric0=eye(n)-G.Y0*G.Y0'/2;
    metrict=eye(n)-G.Yt*G.Yt'/2;
    Gram0=zeros(d);Gramt=zeros(d);
    for k=1:d
        V0=reshape(Z0(:,k),n,p);
        Vt=reshape(Zt(:,k),n,p);
        for q=1:d
            W0=reshape(Z0(:,q),n,p);
            Wt=reshape(Zt(:,q),n,p);
            Gram0(k,q)=trace(V0'*metric0*W0);
            Gramt(k,q)=trace(Vt'*metrict*Wt);
        end
    end
    verifyLessThan(tc,norm(Gram0-eye(d),'fro'),2e-11);
    verifyLessThan(tc,norm(Gramt-eye(d),'fro'),2e-11);
    % Parallel transport of initial velocity must give final velocity.
    K=G.horizontal_geodesic_generator;
    v0=zeros(d,1);
    for q=1:d
        v0(q)=K(G.basis_row(q),G.basis_col(q));
    end
    verifyLessThan(tc,norm(G.U*v0-v0),2e-11);
    Vt=reshape(Zt*(G.U*v0),n,p);
    exact=G.moving_frame*K(:,1:p);
    verifyLessThan(tc,norm(Vt-exact,'fro'),2e-11);
    Bwd=stiefel_parallel_transport(n,p,A,B,-t);
    verifyLessThan(tc,norm(Bwd*G.U-eye(d),'fro'),2e-11);
end
end

function test_sphere_and_full_group_special_cases(tc)
for n=[2 3 6 10]
    B=sin((1:n-1)')*.41;
    Q=stiefel_parallel_transport(n,1,0,B,.9);
    verifyLessThan(tc,norm(Q-eye(n-1),'fro'),1e-12);
end
for n=[2 3 4]
    A=zeros(n);
    A(1,2)=-.41;A(2,1)=.41;
    Q=stiefel_parallel_transport(n,n,A,zeros(0,n),.6);
    verifyLessThan(tc,norm(Q'*Q-eye(n*(n-1)/2),'fro'), ...
        2e-12);
end
end

function test_independent_ambient_levi_civita_ode(tc)
% Independently discretize the canonical-metric Levi-Civita equation,
% NOT the homogeneous connection generator used in the implementation.
n=5;p=2;t=.83;
A=[0 -.32;.32 0];
B=[.27 -.14;.09 .21;-.13 .06];
M=stiefel_parallel_transport(n,p,A,B,t,'model');
K=M.horizontal_geodesic_generator;
options=odeset('RelTol',2e-11,'AbsTol',2e-12);
for idx=[1 3 6]
    F=horizontal_core(n,M.basis_row(idx),M.basis_col(idx));
    V0=F(:,1:p);
    [~,Z]=ode45(@(s,v)ambient_parallel_rhs(s,v,K,n,p), ...
        [0 t],V0(:),options);
    reference=reshape(Z(end,:),n,p);
    E=zeros(n,p);
    for j=1:M.dimension
        Fj=horizontal_core(n,M.basis_row(j),M.basis_col(j));
        E=E+M.U(j,idx)*Fj(:,1:p);
    end
    prediction=M.moving_frame*E;
    verifyLessThan(tc,norm(prediction-reference,'fro'),2e-8);
    verifyLessThan(tc,norm(M.Yt'*reference+reference'*M.Yt, ...
        'fro'),2e-8);
    metric=eye(n)-M.Yt*M.Yt'/2;
    verifyLessThan(tc,abs(trace(reference'*metric*reference)-1), ...
        2e-8);
end
end

function test_spectral_unknown_validation(tc)
S=stiefel_parallel_transport(5,2,[],[],1,'spectral');
verifyEqual(tc,S.dimension,7);
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyEmpty(tc,S.finite.moments);
verifyError(tc,@()stiefel_parallel_transport(1,1), ...
    'stiefel_parallel_transport:Order');
verifyError(tc,@()stiefel_parallel_transport(4,5), ...
    'stiefel_parallel_transport:Order');
verifyError(tc,@()stiefel_parallel_transport(4,2,eye(2)), ...
    'stiefel_parallel_transport:Skew');
verifyError(tc,@()stiefel_parallel_transport(4,2,zeros(2), ...
    ones(3,2)),'stiefel_parallel_transport:Tangent');
verifyError(tc,@()stiefel_parallel_transport(4,2,[],[],1, ...
    'operator'),'stiefel_parallel_transport:Mode');
[U,p]=stiefel_parallel_transport();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(p,'orthogonal')));
end

function F=horizontal_core(n,i,j)
F=zeros(n);
F(i,j)=1;
F(j,i)=-1;
end

function dv=ambient_parallel_rhs(t,v,K,n,p)
% Levi-Civita ODE for g_Y(V,W)=tr(V'*(I-Y*Y'/2)*W).
% Ambient metric tensor is H_Y=I-Y*Y'/2; its inverse is I+Y*Y'.
% Its Christoffel tensor comes from the metric derivative, then the
% tangential projection is removed while differentiating Y'*V+V'*Y=0.
G=expm(t*K);
Y=G(:,1:p);
X=G*K(:,1:p);
V=reshape(v,n,p);
T=Y*(X'*V+V'*X)+X*(Y'*V)+V*(Y'*X) ...
   -X*(V'*Y)-V*(X'*Y);
Gamma=-(T+Y*(Y'*T))/4;
symG=(Y'*Gamma+Gamma'*Y)/2;
symX=(X'*V+V'*X)/2;
dvmat=-Gamma+Y*(symG-symX);
dv=dvmat(:);
end
