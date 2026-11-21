function tests=test_spin_spherical_harmonic_transform
%TEST_SPIN_SPHERICAL_HARMONIC_TRANSFORM Strict block/rectangular geometry.
tests=functiontests(localfunctions);
end

function test_wigner_jacobi_explicit_conventions(tc)
% Independent Wigner matrix exponentials versus normalized Jacobi P.
x=[-.8;-.25;.17;.74];
for ell=0:6
    for m=-ell:ell
        for n=-ell:ell
            G=um_wigner_jacobi(ell+1,m,n,x);
            p=ell-max(abs(m),abs(n))+1;
            for j=1:numel(x)
                D=wigner_rotation(ell,0,acos(x(j)),0);
                d=D(ell+1+m,ell+1+n);
                verifyLessThan(tc,abs(G(j,p)-sqrt((2*ell+1)/2)*d), ...
                    5e-11*(ell+1));
            end
        end
    end
end
end

function test_square_spin_block_gaussian_quadrature(tc)
for L=[1 2 3 5 9 14]
    for s=-min(L-1,2):min(L-1,2)
        for m=-min(L-1,3):min(L-1,3)
            N=L-max(abs(s),abs(m));a=abs(m+s);b=abs(m-s);
            Q=spin_spherical_harmonic_transform(L,s,m);
            verifyEqual(tc,size(Q),[N N]);
            verifyTrue(tc,isreal(Q));
            verifyLessThan(tc,norm(Q'*Q-eye(N),'fro'),8e-11*N);
            reference=(-1)^(s+max(0,m+s))*gauss_jacobi_transform(N,a,b);
            verifyLessThan(tc,norm(Q-reference,'fro'),2e-13*N);
            S=spin_spherical_harmonic_transform(L,s,m,'spectral');
            verifyEqual(tc,S.finite.kind,'unknown');
            verifyEqual(tc,S.dimension,N);
            verifyEqual(tc,S.spin,s);
            P=spin_spherical_harmonic_transform(L,s,m,'parameters');
            verifyEqual(tc,P.fixed_block_dimension,N);
            verifyEqual(tc,P.sampling_rows,L*(2*L-1));
            verifyEqual(tc,P.sampling_columns,L^2-s^2);
            verifyLessThan(tc,abs(sum(P.theta_weights)-2),3e-14*L);
            verifyTrue(tc,all(diff(P.theta_cos_nodes)>0)||L==1);
        end
    end
end
end

function test_rectangular_spin_sampling_isometric(tc)
for L=1:7
    for s=-min(2,L-1):min(2,L-1)
        A=spin_spherical_harmonic_transform(L,s,0,'sampling');
        F=spin_spherical_harmonic_transform(L,s,0,'sampling_operator');
        N=L^2-s^2;M=L*(2*L-1);
        verifyEqual(tc,size(A),[M N]);
        verifyEqual(tc,F([],'size'),[M N]);
        verifyLessThan(tc,norm(A'*A-eye(N),'fro'),7e-11*max(N,1));
        X=randn(N,3)+1i*randn(N,3);
        Y=randn(M,2)+1i*randn(M,2);
        verifyLessThan(tc,norm(F(X,'notransp')-A*X,'fro'), ...
            8e-11*max(N,1));
        verifyLessThan(tc,norm(F(Y,'transp')-A'*Y,'fro'), ...
            8e-11*max(N,1));
        verifyLessThan(tc,norm(F(Y,'ctransp')-A'*Y,'fro'), ...
            8e-11*max(N,1));
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
            8e-11*max(N,1));
        if L>1
            verifyGreaterThan(tc,M,N);
        end
    end
end
end

function test_spin_conjugation_and_scalar_spherical_harmonics(tc)
for L=[2 3 5 7]
    [x,w]=um_legendre_gauss_rule(L);M=2*L-1;
    phi=2*pi*(0:M-1)/M;
    for s=0:min(L-1,2)
        Ap=spin_spherical_harmonic_transform(L,s,0,'sampling');
        Am=spin_spherical_harmonic_transform(L,-s,0,'sampling');
        for ell=s:L-1
            for m=-ell:ell
                ip=ell^2-s^2+m+ell+1;
                im=ell^2-s^2-m+ell+1;
                verifyLessThan(tc,norm(conj(Ap(:,ip))- ...
                    (-1)^(m+s)*Am(:,im)), ...
                    4e-11*L);
            end
        end
    end
    % For s=0,l=0,m=0: the scalar Y_00 weighted sample equals
    % sqrt(wtheta*2pi/M)/sqrt(4pi)=sqrt(wtheta/(2M)).
    A0=spin_spherical_harmonic_transform(L,0,0,'sampling');
    expected=repmat(sqrt(w/(2*M)),M,1);
    verifyLessThan(tc,norm(A0(:,1)-expected),3e-12*L);
    % Independent degree ell=1,m=1 scalar harmonic Condon-Shortley sign.
    if L>=2
        col=1^2+1+1+1;
        sinTheta=sqrt(1-x.^2);
        Y11=-sqrt(3/(8*pi))*sinTheta;
        e=reshape((sqrt(w*2*pi/M).*Y11)*exp(1i*phi),[],1);
        verifyLessThan(tc,norm(A0(:,col)-e),2e-11*L);
    end
end
end

function test_fft_large_and_modes(tc)
L=24;s=2;F=spin_spherical_harmonic_transform(L,s,0,'sampling_operator');
[M,N]=size_contract(F);
verifyEqual(tc,M,L*(2*L-1));
verifyEqual(tc,N,L^2-s^2);
X=randn(N,3)+1i*randn(N,3);
verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
    2e-9*L);
verifyError(tc,@()spin_spherical_harmonic_transform(13,2,0,'sampling'), ...
    'spin_spherical_harmonic_transform:SamplingSize');
verifyError(tc,@()spin_spherical_harmonic_transform(4,4), ...
    'spin_spherical_harmonic_transform:Order');
verifyError(tc,@()spin_spherical_harmonic_transform(4,0,0,'nonsense'), ...
    'spin_spherical_harmonic_transform:Mode');
verifyError(tc,@()F(zeros(N-1,1),'notransp'), ...
    'um_spin_sampling_apply:Dimension');
verifyError(tc,@()F(zeros(M,1),'wrong'), ...
    'um_spin_sampling_apply:Flag');
end

function [M,N]=size_contract(F)
sz=F([],'size');M=sz(1);N=sz(2);
end
