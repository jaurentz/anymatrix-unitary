function tests=test_spherical_mask_slepian
%TEST_SPHERICAL_MASK_SLEPIAN Weighted rotated exact spherical concentration.
tests=functiontests(localfunctions);
end

function test_reduction_to_existing_cap_and_lune_and_empty(tc)
for L=[1 2 4 7]
    n=L^2;
    E=spherical_mask_slepian(L,[]);
    verifyEqual(tc,E,eye(n),'AbsTol',0);
    C=make_component('cap',1,[0 0 0],1.04);
    K=spherical_mask_slepian(L,C,'concentration');
    verifyLessThan(tc,norm(K-spherical_slepian(L,1.04, ...
        'concentration'),'fro'),2e-11*n);
    C=make_component('lune',1,[0 0 0],[.27 1.23]);
    K=spherical_mask_slepian(L,C,'concentration');
    verifyLessThan(tc,norm(K-spherical_lune_slepian( ...
        L,.27,1.23,'concentration'),'fro'), ...
        2e-11*n);
    C=make_component('cap',0,[.3 1.2 .2],1.04);
    verifyEqual(tc,spherical_mask_slepian(L,C),eye(n), ...
        'AbsTol',0);
end
end

function test_arbitrary_3D_rotated_caps_against_independent_integral(tc)
% Independent integration over a rotated spherical cap, evaluating
% real Condon--Shortley harmonics in physical Cartesian coordinates;
% unlike implementation, it never invokes Wigner coefficients,
% spherical_slepian, or an unrotated concentration matrix.
for L=[2 3 5]
    for a=[.27 -.41]
        ang=[a .84 -.27];cap=1.11;
        C=make_component('cap',1,ang,cap);
        K=spherical_mask_slepian(L,C,'concentration');
        Q=direct_rotated_patch(L,'cap',ang,cap);
        verifyLessThan(tc,norm(K-Q,'fro'),4e-10*L^2);
        verifyLessThan(tc,abs(trace(K)-L^2*(1-cos(cap))/2), ...
            5e-10*L^2);
        % Rotation changes eigenvectors, not concentration eigenvalues.
        lambda=eig(K);base=eig(spherical_slepian(L,cap, ...
            'concentration'));
        verifyLessThan(tc,norm(sort(lambda)-sort(base)), ...
            2e-10*L^2);
    end
end
end

function test_rotated_mixed_mask_and_spectrum(tc)
for L=[2 3 5 8]
    n=L^2;
    C=[make_component('cap',.41,[.39 .63 -.27],1.03), ...
       make_component('lune',.34,[-.17 .42 .51],[.22 1.47])];
    P=spherical_mask_slepian(L,C,'parameters');
    K=P.concentration;Q=P.basis;lam=P.concentration_eigenvalues;
    verifyLessThan(tc,norm(K-K','fro'),5e-13*n);
    verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),5e-10*n);
    verifyLessThan(tc,norm(Q'*K*Q-diag(lam),'fro'), ...
        5e-9*n);
    verifyLessThan(tc,abs(trace(K)-P.shannon_number), ...
        6e-10*n);
    verifyLessThanOrEqual(tc,max(lam),1+2e-10);
    verifyGreaterThanOrEqual(tc,min(lam),-2e-10);
    verifyLessThan(tc,norm(spherical_mask_slepian(L,C) ...
        -Q','fro'),2e-12*n);
    verifyLessThan(tc,norm(spherical_mask_slepian(L,C,'basis') ...
        -Q,'fro'),2e-12*n);
    if L<=3
        R=.41*direct_rotated_patch(L,'cap',C(1).angles, ...
            C(1).parameters)+.34*direct_rotated_patch(L, ...
            'lune',C(2).angles,C(2).parameters);
        verifyLessThan(tc,norm(K-R,'fro'),4e-10*n);
    end
    S=spherical_mask_slepian(L,C,'spectral');
    verifyEqual(tc,S.dimension,n);
    verifyEqual(tc,S.finite.kind,'unknown');
    verifyEmpty(tc,S.finite.eigenvalues);
    verifyEmpty(tc,S.finite.moments);
end
end

function test_complementary_caps_and_rotation_conventions(tc)
for L=[2 3 5 8]
    n=L^2;cap=.83;
    north=make_component('cap',1,[0 0 0],cap);
    south=make_component('cap',1,[0 pi 0],pi-cap);
    A=spherical_mask_slepian(L,north,'concentration');
    B=spherical_mask_slepian(L,south,'concentration');
    verifyLessThan(tc,norm(A+B-eye(n),'fro'),8e-10*n);
    C=make_component('cap',1,[.4 .7 -.2],cap);
    C2=make_component('cap',1,[.4 .7 -.2],cap);
    verifyEqual(tc,spherical_mask_slepian(L,C), ...
        spherical_mask_slepian(L,C2),'AbsTol',0);
end
end

function test_validation_and_weighted_not_binary_union(tc)
[U,p]=spherical_mask_slepian();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(p,'orthogonal')));
C=make_component('cap',.6,[0 0 0],1.1);
D=make_component('cap',.5,[0 0 0],1.1);
verifyError(tc,@()spherical_mask_slepian(3,[C,D]), ...
    'spherical_mask_slepian:WeightSum');
C=make_component('cap',.4,[0 0 0],1.1);
D=make_component('cap',.3,[0 0 0],1.1);
K=spherical_mask_slepian(3,[C,D],'concentration');
verifyLessThan(tc,norm(K-.7*spherical_slepian(3,1.1, ...
    'concentration'),'fro'),3e-12);
C=make_component('cap',-.1,[0 0 0],1.1);
verifyError(tc,@()spherical_mask_slepian(3,C), ...
    'spherical_mask_slepian:Weight');
C=make_component('cap',1,[0 0 0],pi+.1);
verifyError(tc,@()spherical_mask_slepian(3,C), ...
    'spherical_mask_slepian:Region');
C=make_component('lune',1,[0 0 0],[.2 2*pi+.01]);
verifyError(tc,@()spherical_mask_slepian(3,C), ...
    'spherical_mask_slepian:Region');
C=make_component('cap',1,[0 0 0],1);
C.angles=[.2 .3];
verifyError(tc,@()spherical_mask_slepian(3,C), ...
    'spherical_mask_slepian:Angles');
C=make_component('other',1,[0 0 0],1);
verifyError(tc,@()spherical_mask_slepian(3,C), ...
    'spherical_mask_slepian:Type');
verifyError(tc,@()spherical_mask_slepian(13), ...
    'MATLAB:notLessEqual');
verifyError(tc,@()spherical_mask_slepian(3,[],'operator'), ...
    'spherical_mask_slepian:Mode');
end

function C=make_component(type,weight,angles,parameters)
C=struct('type',type,'weight',weight, ...
    'angles',angles,'parameters',parameters);
end

function K=direct_rotated_patch(L,type,angles,parameters)
% Physical SO(3) rotation of coordinates, independent of the code's
% real Wigner basis/permutation. Integrate the *actual* real harmonic
% product over the rotated cap or lune using independent quadrature.
n=L^2;
q=3*L+8;
[xx,wx]=gauss_legendre(q);
if strcmp(type,'cap')
    cap=parameters;
    a=cos(cap);b=1;
    x=(a+b)/2+(b-a)*xx/2;
    weightX=wx*(b-a)/2;
    nphi=8*L+9;
    phi=2*pi*(0:nphi-1)'/nphi;
    weightPhi=(2*pi/nphi)*ones(nphi,1);
else
    % Lune products may have sqrt(1-x^2) factors of odd
    % associated-Legendre total order. Integrate smooth colatitude
    % theta with sin(theta) weight, not nonsmooth endpoint x.
    [ztheta,wtheta]=gauss_legendre(4*L+12);
    theta=(pi/2)*(ztheta+1);
    x=cos(theta);weightX=(pi/2)*wtheta.*sin(theta);
    a=parameters(1);b=a+parameters(2);
    [r,w]=gauss_legendre(8*L+12);
    phi=(a+b)/2+(b-a)*r/2;
    weightPhi=(b-a)*w/2;
end
R=physical_rotation(angles);
K=zeros(n);
for i=1:numel(x)
    z=x(i);s=sqrt(max(0,1-z*z));
    for j=1:numel(phi)
        v=R*[s*cos(phi(j));s*sin(phi(j));z];
        theta=acos(max(-1,min(1,v(3))));
        az=atan2(v(2),v(1));
        Y=real_harmonics(L,theta,az);
        K=K+weightX(i)*weightPhi(j)*(Y*Y');
    end
end
end

function Y=real_harmonics(L,theta,az)
Y=zeros(L^2,1);
idx=0;
for m=0:L-1
    if m==0,kinds=0;else,kinds=[1 2];end
    for kind=kinds
        for ell=m:L-1
            P=legendre(ell,cos(theta));
            b=sqrt((2*ell+1)/2*exp(gammaln(ell-m+1) ...
                -gammaln(ell+m+1)))*P(m+1);
            if kind==0,v=1/sqrt(2*pi);
            elseif kind==1,v=cos(m*az)/sqrt(pi);
            else,v=sin(m*az)/sqrt(pi);end
            idx=idx+1;Y(idx)=b*v;
        end
    end
end
end

function R=physical_rotation(angles)
a=angles(1);b=angles(2);g=angles(3);
Za=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
Yb=[cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b)];
Zg=[cos(g),-sin(g),0;sin(g),cos(g),0;0,0,1];
R=Za*Yb*Zg;
end

function [x,w]=gauss_legendre(n)
j=(1:n-1)';
off=j./sqrt(4*j.^2-1);
J=diag(off,1)+diag(off,-1);
[V,D]=eig(J);
[x,ind]=sort(real(diag(D)));
w=2*real(V(1,ind).').^2;
end
