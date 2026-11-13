%UNIT_CIRCLE_SPECTRAL_BANDS Finite eigenphases, gaps, and limiting spectral bands.
% Run from the repository root after: addpath(pwd)
% Set plotSpectralBands=false before running to skip the optional figure.
% MATLAB R2020b or later, no additional toolboxes.
%
% IMPORTANT: a finite n-by-n unitary has n eigenvalues (an atomic
% empirical trace measure). Analytic continuous densities and CMV bands
% below are LIMITING/INFINITE-OPERATOR references, not finite densities.
if ~exist('plotSpectralBands','var') || isempty(plotSpectralBands)
    plotSpectralBands=true;
end

% 1. TWO PRESCRIBED ARCS, exact finite eigenphases, dense mixing.
% This tests a hard spectral gap without relying on an asymptotic model.
phaseBands=[linspace(-2.70,-1.30,34),linspace(0.35,1.35,22)]';
nPrescribed=numel(phaseBands);
U=phase_circulant(nPrescribed,phaseBands);
Uop=phase_circulant(nPrescribed,phaseBands,'operator');
v=exp(1i*(0:nPrescribed-1)'*.23);
assert(norm(U'*U-eye(nPrescribed),'fro')<1e-10)
assert(norm(Uop(v,'notransp')-U*v)<1e-10)
for r=1:3
    assert(abs(trace(U^r)-sum(exp(1i*r*phaseBands)))<1e-9)
end

% 2. DIRAC AUTOMATON: two analytic bands, known limiting density.
nDirac=96;
mass=0.45;
D=dirac_automaton(nDirac,mass);
nu=sqrt(1-mass^2);
k=(0:nDirac-1)';
omegaD=acos(nu*cos(2*pi*k/nDirac));
thetaD=angle(eig(full(D)));
edgeD=asin(mass);
% Allowed bands on [-pi,pi]:
%   [edgeD,pi-edgeD] and [-pi+edgeD,-edgeD].
assert(all(abs(thetaD)>=edgeD-2e-12))
assert(all(abs(thetaD)<=pi-edgeD+2e-12))
assert(norm(full(D'*D-speye(2*nDirac)),'fro')<1e-10)
for r=1:3
    assert(abs(trace(D^r)-2*sum(cos(r*omegaD)))<1e-9)
end
% Limiting equal-weight eigenphase density:
% rho(theta)=abs(sin(theta))/(2*pi*sqrt(nu^2-cos(theta)^2))
% on abs(cos(theta))<nu, zero elsewhere. Both bands carry total mass 1.
thetaGrid=linspace(-pi,pi,1601)';
rhoD=zeros(size(thetaGrid));
insideD=abs(cos(thetaGrid))<nu;
rhoD(insideD)=abs(sin(thetaGrid(insideD)))./ ...
    (2*pi*sqrt(nu^2-cos(thetaGrid(insideD)).^2));

% 3. PERIOD-2 CMV: explicit infinite-operator bands via the discriminant.
% For real Verblunsky values a0,a1 and their rho_j=sqrt(1-a_j^2),
% Delta(e^{i theta})=2*(cos(theta)+a0*a1)/(rho0*rho1).
% Infinite bands are exactly |Delta|<=2. Finite CMV eigenvalues can
% include boundary-condition outliers in these infinite-operator gaps.
a0=0.48;a1=-0.31;
pattern=[a0;a1];
rhoProduct=sqrt((1-a0^2)*(1-a1^2));
discP=2*(cos(thetaGrid)+a0*a1)/rhoProduct;
inBandP=abs(discP)<=2;
cosLower=-rhoProduct-a0*a1;
cosUpper= rhoProduct-a0*a1;
assert(cosLower>-1 && cosUpper<1)
positiveEdgesP=sort(acos([cosLower,cosUpper]));
% Two infinite bands: theta in +/-[positiveEdgesP(1),positiveEdgesP(2)].
P=periodic_cmv(96,pattern);
thetaP=angle(eig(full(P)));
Pop=periodic_cmv(96,pattern,'cmv',1,[],'operator');
vP=exp(1i*.12*(0:95)');
assert(norm(Pop(vP,'notransp')-P*vP)<1e-10)
numBoundaryLevels=nnz(abs(2*(cos(thetaP)+a0*a1)/rhoProduct)>2+1e-9);
% A finite boundary level in a gap does not contradict the bulk bands.

% 4. GERONIMUS CMV: constant real Verblunsky coefficient.
% As a special period-two case a0=a1=a, the infinite essential support
% is |theta|>=2*asin(abs(a)), a proper unit-circle arc.
% Depending on the measure/closure there may also be a point mass or
% finite boundary eigenvalues inside the essential-spectrum gap.
a=0.36;
G=geronimus_cmv(96,a);
thetaG=angle(eig(full(G)));
geronimusEdge=2*asin(abs(a));
assert(norm(full(G'*G-speye(96)),'fro')<1e-10)

% 5. CYCLE PROPAGATOR: a wrapped cosine dispersion with exact phase law.
% Choose 2*t<pi so the limiting arcsine density has one circular arc
% [-2*t,2*t] with no phase wrapping.
tCycle=0.65;
nCycle=96;
T=cycle_propagator(nCycle,tCycle);
thetaT=angle(eig(T));
momentum=2*pi*(0:nCycle-1)'/nCycle;
exactPhases=-2*tCycle*cos(momentum);
assert(all(abs(thetaT)<=2*tCycle+1e-11))
for r=1:3
    assert(abs(trace(T^r)-sum(exp(1i*r*exactPhases)))<1e-9)
end
rhoT=zeros(size(thetaGrid));
insideT=abs(thetaGrid)<2*tCycle;
rhoT(insideT)=1./(pi*sqrt((2*tCycle)^2-thetaGrid(insideT).^2));

% 6. SPLIT-STEP WALK: two exact quasienergy branches (uniform periodic
% coins); the finite spectrum is the momentum sample of these branches.
nWalk=72;
theta1=0.3;theta2=0.7;
S=split_step_walk(nWalk,theta1,theta2);
thetaS=angle(eig(full(S)));
kw=2*pi*(0:nWalk-1)'/nWalk;
omegaS=acos(cos(theta1)*cos(theta2)*cos(kw)-sin(theta1)*sin(theta2));
for r=1:3
    assert(abs(trace(S^r)-2*sum(cos(r*omegaS)))<1e-9)
end

% 7. Circular ensembles, unitary dilations, and OPUC quadrature.
% CUE beta=2 is an ensemble law, not a deterministic finite density.
C=circular_beta(64,2,42);
assert(norm(full(C'*C-speye(64)),'fro')<1e-11)

% Unitary dilation of a strict contraction.
contraction=.3*full(cycle_shift(8));
J=julia_unitary(8,contraction);
assert(norm(J'*J-eye(16),'fro')<1e-11)

% Poisson/Szego continuous CYCLIC-VECTOR weight -> finite CMV quadrature.
% The Gauss-Szego weights are not the equal-weight trace spectral measure.
mu=szego_bernstein_measure(.3+.2i);
alpha=szego_verblunsky(mu,19);
B=szego_matrix(mu,20,'cmv',1);
[z,w]=szego_quadrature(B);
assert(numel(alpha)==19 && abs(sum(w.*z)-conj(mu.alpha0))<1e-10)

fprintf('Dirac band edges +/-[%.4f, %.4f] rad.\n',edgeD,pi-edgeD);
fprintf('Period-2 CMV infinite bands +/-[%.4f, %.4f] rad; %d finite gap levels.\n', ...
        positiveEdgesP(1),positiveEdgesP(2),numBoundaryLevels);
fprintf('Geronimus essential-spectrum gap |theta|<%.4f rad.\n',geronimusEdge);
fprintf('Tight-binding limiting arc [-%.4f, %.4f] rad.\n',2*tCycle,2*tCycle);
fprintf('Processed %d finite eigenphases across six band examples.\n', ...
 nPrescribed+2*nDirac+numel(thetaP)+numel(thetaG)+nCycle+2*nWalk);

% Visualization distinguishes finite eigenphases from analytic support.
if plotSpectralBands
    figure('Name','Unitary spectral bands and gaps');
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    nexttile
    scatter(phaseBands,zeros(size(phaseBands)),12,'filled');
    xlim([-pi pi]);ylim([-.1 .3]);title('Prescribed finite arcs');
    xlabel('Eigenphase');ylabel('Finite eigenvalues');

    nexttile
    histogram(thetaD,40,'Normalization','pdf');
    hold on
    plot(thetaGrid,rhoD,'LineWidth',1.5);
    xlim([-pi pi]);title('Dirac: exact limiting density');
    xlabel('Eigenphase');ylabel('Density');
    legend('Finite eigenphases','Limiting density');

    nexttile
    scatter(thetaP,.15*ones(size(thetaP)),10,'filled');
    hold on
    plot(thetaGrid,double(inBandP),'LineWidth',1.4);
    xlim([-pi pi]);ylim([-.12 1.2]);
    title('Period-2 CMV: infinite bands');
    xlabel('Eigenphase');ylabel('Allowed band = 1');
    legend('Finite eigenphases','Infinite band indicator');

    nexttile
    scatter(thetaG,.15*ones(size(thetaG)),10,'filled');
    hold on
    plot(thetaGrid,double(abs(thetaGrid)>=geronimusEdge),'LineWidth',1.4);
    xlim([-pi pi]);ylim([-.12 1.2]);
    title('Geronimus: essential arc and gap');
    xlabel('Eigenphase');ylabel('Essential support = 1');

    nexttile
    histogram(thetaT,35,'Normalization','pdf');
    hold on
    plot(thetaGrid,rhoT,'LineWidth',1.5);
    xlim([-pi pi]);title('Tight binding: arcsine limit');
    xlabel('Eigenphase');ylabel('Density');

    nexttile
    scatter(thetaS,.1*ones(size(thetaS)),10,'filled');
    hold on
    scatter([-omegaS;omegaS],.3*ones(2*nWalk,1),8);
    xlim([-pi pi]);ylim([-.1 .5]);
    title('Split-step walk: exact branches');
    xlabel('Eigenphase');ylabel('Finite / exact sample');
end

