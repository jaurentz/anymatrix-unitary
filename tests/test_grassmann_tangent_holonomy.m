function tests=test_grassmann_tangent_holonomy
%TEST_GRASSMANN_TANGENT_HOLONOMY Canonical Grassmann fiber transport.
tests=functiontests(localfunctions);
end

function test_intrinsic_dimensions_and_isometry_real_complex(tc)
for field={'real','complex'}
    mode=field{1};
    for nr={[3 1],[4 2],[5 2],[6 3]}
        sizes=nr{1};n=sizes(1);r=sizes(2);
        frames=make_frames(n,r,mode);
        U=grassmann_tangent_holonomy(n,r,frames,mode);
        D=grassmann_tangent_holonomy(n,r,frames,mode,'model');
        dim=r*(n-r);
        verifyEqual(tc,size(U),[dim dim]);
        verifyEqual(tc,D.tangent_dimension,dim);
        verifyLessThan(tc,norm(U'*U-eye(dim),'fro'),5e-10*dim);
        verifyLessThan(tc,norm(D.U-U,'fro'),5e-11*dim);
        verifyLessThan(tc,D.unitarity_residual,5e-10*dim);
        if strcmp(mode,'real')
            verifyTrue(tc,isreal(U));
            verifyLessThan(tc,norm(U.'*U-eye(dim),'fro'), ...
                5e-10*dim);
        end
        for k=1:size(frames,3)
            Y=frames(:,:,k);
            A=D.ambient_lifts(:,:,k);
            verifyLessThan(tc,norm(A'*A-eye(n),'fro'),5e-10*n);
            verifyLessThan(tc,norm(Y'*D.geodesic_initial_velocities(:,:,k),'fro'), ...
                5e-10*n);
            verifyLessThan(tc,norm(D.links(:,:,k)'*D.links(:,:,k) ...
                   -eye(dim),'fro'),5e-10*dim);
            verifyLessThan(tc,norm(D.alignments(:,:,k)' ...
                   *D.alignments(:,:,k)-eye(r),'fro'),5e-10);
        end
    end
end
end

function test_inverse_geodesic_polygon_and_gauge_invariance(tc)
for field={'real','complex'}
    mode=field{1};n=5;r=2;
    frames=make_frames(n,r,mode);
    U=grassmann_tangent_holonomy(n,r,frames,mode);
    inverse=frames(:,:,[1 3 2]);
    V=grassmann_tangent_holonomy(n,r,inverse,mode);
    dim=r*(n-r);
    verifyLessThan(tc,norm(U*V-eye(dim),'fro'),8e-10*dim);
    % Arbitrary unitary gauge at each vertex changes the input frame
    % but leaves intrinsic tangent holonomy eigenphases invariant.
    G=frames;
    for j=1:3
        A=zeros(r);
        if r>=2,A(1,2)=.13*j;A(2,1)=-.13*j;end
        Q=expm(A);
        if strcmp(mode,'complex'),Q=exp(.19i*j)*Q;end
        G(:,:,j)=frames(:,:,j)*Q;
    end
    W=grassmann_tangent_holonomy(n,r,G,mode);
    ev=eig(U);ew=eig(W);
    verifyLessThan(tc,max(arrayfun(@(z)min(abs(ev-z)),ew)),3e-8);
    verifyLessThan(tc,max(arrayfun(@(z)min(abs(ew-z)),ev)),3e-8);
end
end

function test_CP1_Fubini_Study_curvature_and_Berry_square(tc)
n=2;r=1;
for phase=[.25 .7 1.2]
    F=make_cp_polygon([1],phase);
    B=grassmann_holonomy(2,1,F);
    U=grassmann_tangent_holonomy(n,r,F,'complex');
    S=grassmann_tangent_holonomy(n,r,F,'complex','spectral');
    verifyEqual(tc,size(U),[1 1]);
    verifyLessThan(tc,abs(U-conj(B)^2),6e-11);
    verifyEqual(tc,S.finite.kind,'exact');
    verifyLessThan(tc,abs(S.finite.eigenvalues()-U),6e-11);
    for k=-3:4
        verifyLessThan(tc,abs(S.finite.moments(k)-U^k),6e-11);
    end
    % With standard CP1 Fubini-Study normalization (curvature 4),
    % the Berry line phase is +/-2*FS area and tangent +/-4*FS
    % area. These comparisons are modulo 2*pi.
    fsArea=abs(angle(B))/2;
    verifyLessThan(tc,abs(abs(angle(U))-4*fsArea),3e-10);
end
end

function test_CP2_genuine_U2_nonabelian_tangent_holonomy(tc)
d1=[1;0];d2=[1;1]/sqrt(2);
F1=make_cp_polygon(d1,.9);
F2=make_cp_polygon(d2,.9);
A=grassmann_tangent_holonomy(3,1,F1,'complex');
B=grassmann_tangent_holonomy(3,1,F2,'complex');
verifyEqual(tc,size(A),[2 2]);
verifyLessThan(tc,norm(A'*A-eye(2),'fro'),3e-11);
verifyLessThan(tc,norm(B'*B-eye(2),'fro'),3e-11);
verifyGreaterThan(tc,norm(A*B-B*A,'fro'),5e-4);
S=grassmann_tangent_holonomy(3,1,F1,'complex','spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
end

function test_independent_connection_ODE_and_parallel_velocity(tc)
for field={'real','complex'}
    mode=field{1};n=5;r=2;
    F=make_frames(n,r,mode);
    for k=1:3
        Y=F(:,:,k);
        Z=F(:,:,mod(k,3)+1);
        [T,A,Ralign,S,Delta,Q0,Q1]= ...
            um_grassmann_tangent_leg(Y,Z);
        Omega=Delta*Y'-Y*Delta';
        verifyLessThan(tc,norm(Omega+Omega','fro'),4e-11);
        verifyLessThan(tc,norm(A*Y-Z*Ralign,'fro'),2e-9);
        verifyLessThan(tc,norm(Q1'*A*Q0-S,'fro'),2e-11);
        verifyLessThan(tc,norm(T-kron(conj(Ralign),S),'fro'), ...
            2e-11);
        Bs=reshape((1:(n-r)*r),n-r,r)/20;
        if strcmp(mode,'complex'),Bs=Bs+1i*Bs/3;end
        V=Q0*Bs;
        reference=parallel_ODE_RK4(Y,Omega,V,180);
        verifyLessThan(tc,norm(reference-A*V,'fro'),5e-8);
        for t=[0 .3 .8 1]
            At=expm(t*Omega);Yt=At*Y;Vt=At*V;
            dot=Omega*Vt;
            verifyLessThan(tc,norm((eye(n)-Yt*Yt')*dot,'fro'), ...
                1e-10);
            verifyLessThan(tc,norm(Yt'*Vt,'fro'),1e-10);
            verifyLessThan(tc,norm(Vt'*Vt-V'*V,'fro'),1e-10);
        end
        % Along the geodesic, its own initial velocity Delta
        % is Levi-Civita parallel.
        endVel=A*Delta;
        verifyLessThan(tc,norm((eye(n)-Z*Z')*Omega*endVel,'fro'), ...
            1e-10);
    end
end
end

function test_validation_and_inventory(tc)
[A,props]=grassmann_tangent_holonomy();
verifyEmpty(tc,A);
verifyTrue(tc,any(strcmp(props,'unitary')));
verifyTrue(tc,any(strcmp(props,'normal')));
verifyEqual(tc,size(grassmann_tangent_holonomy(4,2)),[4 4]);
verifyError(tc,@()grassmann_tangent_holonomy(1), ...
    'grassmann_tangent_holonomy:Order');
verifyError(tc,@()grassmann_tangent_holonomy(3,3), ...
    'grassmann_tangent_holonomy:Rank');
F=make_frames(3,1,'complex');
verifyError(tc,@()grassmann_tangent_holonomy(3,1,F,'real'), ...
    'grassmann_tangent_holonomy:Frames');
verifyError(tc,@()grassmann_tangent_holonomy(3,1,F,'complex','foo'), ...
    'grassmann_tangent_holonomy:Mode');
F(:,:,2)=F(:,:,1);
verifyError(tc,@()grassmann_tangent_holonomy(3,1,F), ...
    'grassmann_tangent_holonomy:CutLocus');
end

function F=make_frames(n,r,field)
F=zeros(n,r,3);
F(:,:,1)=eye(n,r);
F(:,:,2)=eye(n,r);
F(1,1,2)=cos(.38);F(r+1,1,2)=sin(.38);
F(:,:,3)=eye(n,r);
if r==1,j=1;else,j=2;end
F(j,j,3)=cos(.53);
if strcmp(field,'complex')
    F(r+1,j,3)=exp(.7i)*sin(.53);
else
    F(r+1,j,3)=sin(.53);
end
end

function F=make_cp_polygon(direction,phase)
direction=direction(:)/norm(direction);
n=numel(direction)+1;
F=zeros(n,1,3);
F(:,1,1)=[1;zeros(n-1,1)];
v=[0;direction];
F(:,1,2)=cos(.5)*F(:,1,1)+sin(.5)*v;
F(:,1,3)=cos(.6)*F(:,1,1)+exp(1i*phase)*sin(.6)*v;
end

function V=parallel_ODE_RK4(Y,Omega,V,steps)
% Horizontal connection equation V'=-Y(t)*(Ydot(t)'*V),
% Y(t)=exp(t Omega)*Y, distinct from the closed-form ambient lift.
h=1/steps;
for j=1:steps
    t=(j-1)*h;
    k1=rhs(Y,Omega,t,V);
    k2=rhs(Y,Omega,t+h/2,V+h*k1/2);
    k3=rhs(Y,Omega,t+h/2,V+h*k2/2);
    k4=rhs(Y,Omega,t+h,V+h*k3);
    V=V+h*(k1+2*k2+2*k3+k4)/6;
end
end

function dV=rhs(Y,Omega,t,V)
Yt=expm(t*Omega)*Y;
dot=Omega*Yt;
dV=-Yt*(dot'*V);
end
