%VIBRON_ESQPT_SECTOR_GAPS Angular-momentum clustering on the unit circle.
% Based on Khalouf-Rivera et al., Phys. Rev. A 109 (2024), 062219.
% Run from repository root after addpath(pwd). No toolbox required.
%
% The U(3) two-dimensional vibron model has one tridiagonal Hamiltonian
% for every conserved angular momentum ell. The difference between the
% low-lying ell=0 and ell=1 energies decays algebraically with N in
% the broken-symmetry regime. The U(2) Lipkin parity-sector gap can
% decay exponentially, often reaching floating-point resolution.
% Do not interpret wrapped finite eigenphases as a continuous density.

%% A physical unitary matrix and its two angular-momentum sectors
N=60;
xi=.5;                    % ESQPT regime: xi > 1/5
time=.12;                 % Small enough to resolve low-energy phase pairs
H0=vibron_2d_esqpt(N,xi,time,0,'hamiltonian');
H1=vibron_2d_esqpt(N,xi,time,1,'hamiltonian');
U0=vibron_2d_esqpt(N,xi,time,0);
U1=vibron_2d_esqpt(N,xi,time,1);
assert(norm(U0'*U0-eye(size(U0,1)),'fro')<1e-10);
assert(norm(U1'*U1-eye(size(U1,1)),'fro')<1e-10);

E0=sort(real(eig(full(H0))));
E1=sort(real(eig(full(H1))));
gap=abs(E0(1)-E1(1));
fprintf('2D vibron, N=%d, xi=%g, ell=0 versus 1: gap=%.4e.\n', ...
        N,xi,gap);

%% Finite eigenphases of the two sectors
% Each eigenvalue lies on the unit circle. Compare matched low-lying
% energy levels, rather than numerically sorting wrapped eigenphases.
ph0=exp(-1i*time*E0);
ph1=exp(-1i*time*E1);
figure;
plot(real(ph0),imag(ph0),'o','MarkerSize',5);
hold on
plot(real(ph1),imag(ph1),'x','MarkerSize',5);
theta=linspace(-pi,pi,400);
plot(cos(theta),sin(theta),':');
axis equal
xlabel('Real part');ylabel('Imaginary part');
title(sprintf('2D vibron: N=%d, xi=%.2f, t=%.2f',N,xi,time));
legend('ell=0','ell=1','unit circle','Location','best');

%% Compare gap decay: U(3) vibron versus U(2) Lipkin
Ns=[12 16 20 24 30 40 60 80 100];
vibronGap=zeros(size(Ns));lipkinGap=zeros(size(Ns));
for j=1:numel(Ns)
    n=Ns(j);
    h0=vibron_2d_esqpt(n,xi,time,0,'hamiltonian');
    h1=vibron_2d_esqpt(n,xi,time,1,'hamiltonian');
    e0=eig(full(h0));
    e1=eig(full(h1));
    vibronGap(j)=abs(min(real(e0))-min(real(e1)));

    % At alpha=0, the existing extended Lipkin generator agrees with
    % the simpler U(2) LMG Hamiltonian in the paper up to a scalar shift.
    H=extended_lipkin_propagator(n,xi,0,time,'hamiltonian');
    even=eig(full(H(1:2:end,1:2:end)));
    odd=eig(full(H(2:2:end,2:2:end)));
    lipkinGap(j)=abs(min(real(even))-min(real(odd)));
end

figure;
semilogy(Ns,vibronGap,'o-','LineWidth',1.2);hold on
semilogy(Ns,max(lipkinGap,eps),'s--','LineWidth',1.2);
xlabel('Number of bosons N');
ylabel('Lowest cross-sector energy gap');
legend('U(3): ell=0 versus ell=1', ...
       'U(2): even versus odd parity','Location','best');
title('ESQPT spectral clustering: algebraic versus exponential');
grid on
% The Lipkin curve saturates at floating-point resolution, not a
% physically meaningful nonzero limiting gap.
