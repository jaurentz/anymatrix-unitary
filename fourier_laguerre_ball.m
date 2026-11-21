function [U,properties]=fourier_laguerre_ball(P,L,tau,output)
%FOURIER_LAGUERRE_BALL Exact square axisymmetric Fourier-Laguerre block.
%   U=FOURIER_LAGUERRE_BALL(P,L,TAU,OUTPUT) is a real orthogonal
%   P*L-by-P*L transform on the AXISYMMETRIC (m=0) bandlimited
%   Fourier--Laguerre subspace of the ball. Radial orders p=0,...,P-1,
%   angular degrees ell=0,...,L-1. P,L positive integers, TAU>0
%   physical radial scale (default 1).
%
%   Radial orthonormal basis in L2(r^2 dr) on [0,infinity):
%      K_p(r)=sqrt(p!/(p+2)!)*TAU^(-3/2)*exp(-r/(2*TAU))
%             *L_p^2(r/TAU).
%   With P generalized Gauss--Laguerre (alpha=2) nodes x_j, ordinary
%   quadrature weights w_j for x^2 exp(-x), sample radii r_j=TAU*x_j
%   and PHYSICAL radial weights v_j=TAU^3*exp(x_j)*w_j,
%      Qr(j,p+1)=sqrt(v_j)*K_p(r_j).
%   The associated positive-leading Gauss-Laguerre transform has
%   sign (-1)^p relative to the standard L_p^2 convention.
%
%   Axisymmetric real Y_ell,0 = sqrt((2ell+1)/(4*pi))*P_ell(cos theta),
%   with L-point Gauss--Legendre quadrature nodes x=cos theta,
%   weights wtheta and analytic phi integration, yields
%      Qa(j,ell+1)=sqrt(2*pi*wtheta_j)*Y_ell,0(theta_j).
%   Then U=kron(Qa,Qr) is exactly orthogonal by product quadrature.
%   Coefficients use radial index fastest. TAU changes the PHYSICAL
%   sampling radii, not the normalized finite orthogonal matrix or
%   its unknown eigenphases.
%
%   'matrix' (default): explicit square orthogonal transform.
%   'operator': Kronecker modewise action without materializing U
%               (O(P*L*(P+L)) per RHS after O(P^3+L^3) setup).
%   'parameters': returns Qr,Qa, nodes, physical weights and sizes.
%   'spectral': unknown finite trace eigenphase law; Gauss nodes
%               are NOT eigenvalues of U.
%   'sampling': FULL-SPHERE rectangular weighted analysis/synthesis
%               E=kron(Aangular,Qr), size [P*L*(2L-1),P*L^2],
%               using L Gauss--Legendre colatitudes and 2L-1
%               equispaced azimuths. E'*E=I, but E is NOT square
%               unitary except L=1. This is a separate isometric
%               sampling interface, not a new matrix-family count.
%               Restricted to P*L^2<=512 to avoid dense allocations.
%
%   This is not a spherical Gauss--Laguerre (Gaussian r^2) transform,
%   nor an assertion that general 3D sampling is square unitary.
%
%   [157] Leistedt and McEwen, IEEE TSP 60 (2012).
%   https://doi.org/10.1109/TSP.2012.2215030
%   [158] McEwen and Leistedt, SampTA (2013).
%   https://arxiv.org/abs/1307.1307
properties={'square','unitary','orthogonal','real','normal', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(P,{'numeric'},{'scalar','integer','positive','finite','<=',128});
if nargin<2||isempty(L),L=3;end
validateattributes(L,{'numeric'},{'scalar','integer','positive','finite','<=',128});
if nargin<3||isempty(tau),tau=1;end
validateattributes(tau,{'numeric'},{'scalar','real','positive','finite'});
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('fourier_laguerre_ball:Mode','OUTPUT must be scalar text.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','operator','parameters','spectral','sampling'}))
    error('fourier_laguerre_ball:Mode','Unknown output mode.');
end
N=double(P)*double(L);
if strcmp(output,'spectral')
    U=unitary_spectral_reference('fourier_laguerre_ball',N);
    U.radial_bandlimit=P;U.angular_bandlimit=L;U.radial_scale=tau;
    U.reference='Leistedt--McEwen, IEEE TSP 60 (2012)';
    U.notes=['Exact square axisymmetric Fourier-Laguerre quadrature ' ...
        'transform; its eigenphases are unknown. Laguerre roots and ' ...
        'Legendre roots are quadrature nodes, NOT eigenphases. ' ...
        'Full 3D angular sampling is generally rectangular isometric.'];
    return
end
Qr=gauss_laguerre_transform(P,2);
Qr=Qr*diag((-1).^(0:P-1));
Qa=gauss_jacobi_transform(L,0,0);
if strcmp(output,'operator')
    U=@(X,flag)um_tensor_apply(X,flag,{Qa,Qr});
    return
end
if strcmp(output,'matrix')
    if N>2048
        error('fourier_laguerre_ball:Size', ...
            'Explicit square matrix is capped at order 2048; use operator.');
    end
    U=kron(Qa,Qr);return
end
% Gauss-Laguerre J for x^2*exp(-x) and Gauss-Legendre J for dx.
pr=(0:P-1)'; jr=(1:P-1)';
Jr=diag(2*pr+3)+diag(sqrt(jr.*(jr+2)),1) ...
    +diag(sqrt(jr.*(jr+2)),-1);
[Vr,Dr]=eig(Jr);
[xr,ir]=sort(diag(Dr));wr=2*(Vr(1,ir).').^2;
% Raw positive physical weights may overflow at very large P or tau.
% They are included only in diagnostic output, never used to form Qr.
rw=tau^3*exp(xr).*wr;
ja=(1:L-1)';ea=ja./sqrt(4*ja.^2-1);
Ja=diag(ea,1)+diag(ea,-1);
[Va,Da]=eig(Ja);
[xa,ia]=sort(diag(Da));wa=2*(Va(1,ia).').^2;
if strcmp(output,'parameters')
    U=struct('radial_transform',Qr,'angular_transform',Qa, ...
        'radial_nodes',tau*xr,'radial_quadrature_weights',rw, ...
        'laguerre_nodes',xr,'laguerre_weighted_weights',wr, ...
        'cos_theta_nodes',xa,'theta_weights',wa, ...
        'radial_bandlimit',P,'angular_bandlimit',L, ...
        'radial_scale',tau,'dimension',N, ...
        'full_spherical_sampling_rows',P*L*(2*L-1), ...
        'full_spherical_coefficients',P*L^2, ...
        'angular_zonal_only',true);
    return
end
if P*L^2>512
    error('fourier_laguerre_ball:SamplingSize', ...
        'Rectangular sampling restricted to P*L^2<=512.');
end
Nphi=2*L-1;
phi=2*pi*(0:Nphi-1)/Nphi;
Aang=zeros(L*Nphi,L^2);
col=0;
for ell=0:L-1
    Plm=legendre(ell,xa.');
    for m=-ell:ell
        col=col+1;
        mm=abs(m);
        C=sqrt((2*ell+1)/(4*pi)* ...
            exp(gammaln(ell-mm+1)-gammaln(ell+mm+1)));
        radial=Plm(mm+1,:).';
        if m==0
            angular=C*ones(1,Nphi);
        elseif m>0
            angular=sqrt(2)*C*cos(mm*phi);
        else
            angular=sqrt(2)*C*sin(mm*phi);
        end
        vals=radial*angular;
        weighted=(sqrt(wa)*ones(1,Nphi)).*vals*sqrt(2*pi/Nphi);
        Aang(:,col)=weighted(:);
    end
end
U=kron(Aang,Qr); % rectangular isometry, explicitly NOT unitary
end
