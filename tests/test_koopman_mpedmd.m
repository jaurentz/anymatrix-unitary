function tests=test_koopman_mpedmd
%TEST_KOOPMAN_MPEDMD Measured Gram-Procrustes Koopman unitary models.
tests=functiontests(localfunctions);
end

function test_default_exact_circle_rotation_spectrum(tc)
for n=[1 2 3 5 11 21]
    Q=koopman_mpedmd(n);
    P=koopman_mpedmd(n,[],[],[],'metric');
    data=koopman_mpedmd(n,[],[],[],'data');
    S=koopman_mpedmd(n,[],[],[],'spectral');
    exact=exp(1i*.37*(0:n-1)');
    verifyEqual(tc,size(Q),[n n]);
    verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),1e-10*n);
    verifyLessThan(tc,norm(Q-diag(exact),'fro'),2e-10*n);
    verifyEqual(tc,S.finite.kind,'exact');
    verifyLessThan(tc,norm(S.finite.eigenvalues()-exact), ...
        1e-14*n);
    for k=[-5 -1 0 1 3 7]
        verifyLessThan(tc,abs(S.finite.moments(k) ...
            -mean(exact.^k)),1e-12);
        verifyLessThan(tc,abs(trace(Q^abs(k))/n ...
            -S.finite.moments(abs(k))),2e-10*n);
    end
    verifyLessThan(tc,norm(P.G-eye(n),'fro'),2e-10*n);
    verifyLessThan(tc,norm(P.K-Q,'fro'),2e-10*n);
    verifyEqual(tc,P.cross_rank,n);
    verifyFalse(tc,P.rank_deficient_completion);
    verifyEqual(tc,size(data.PX,2),n);
    verifyEqual(tc,size(data.PY),size(data.PX));
    verifyEqual(tc,numel(data.snapshot_x),size(data.PX,1));
    verifyLessThan(tc,abs(sum(data.weights)-1),2e-13);
end
end

function test_author_mpEDMD_Gram_SVD_orientation(tc)
% Directly reimplement Colbrook's original mpEDMD.m orientation:
% SVD of G^(-1/2)*A'*G^(-1/2) gives U,V, and Q=V*U'.
n=4;p=47;
x=(2*pi/p)*(0:p-1)';
y=mod(x+.37+.13*sin(3*x),2*pi);
PX=[ones(p,1),cos(x),sin(x),cos(2*x)];
PY=[ones(p,1),cos(y),sin(y),cos(2*y)];
w=1+.37*sin(5*x).^2;
data=koopman_mpedmd(n,PX,PY,w,'metric');
W=w/sum(w);
G=PX'*(W.*PX);
A=PX'*(W.*PY);
G=(G+G')/2;
[V,D]=eig(G);
Ghalf=V*diag(sqrt(diag(D)))*V';
Ginv=V*diag(1./sqrt(diag(D)))*V';
[Us,~,Vs]=svd(Ginv*A'*Ginv);
Qauthor=Vs*Us';
Kauthor=Ginv*Qauthor*Ghalf;
verifyLessThan(tc,norm(data.G-G,'fro'),2e-12);
verifyLessThan(tc,norm(data.A-A,'fro'),2e-12);
verifyLessThan(tc,norm(data.Q-Qauthor,'fro'),5e-10);
verifyLessThan(tc,norm(data.K-Kauthor,'fro'),5e-10);
verifyLessThan(tc,norm(data.H'*data.H-G,'fro'),5e-11);
verifyLessThan(tc,norm(data.Q'*data.Q-eye(n),'fro'),5e-11);
verifyLessThan(tc,norm(data.K'*G*data.K-G,'fro'),5e-10);
verifyLessThan(tc,norm(data.EDMD-G\A,'fro'),5e-11);
verifyGreaterThan(tc,norm(data.EDMD'*G*data.EDMD-G,'fro'),1e-3);
verifyGreaterThan(tc,norm(data.K'*data.K-eye(n),'fro'),1e-3);
end

function test_noninvariant_noisy_dictionary_and_analytic_unknown(tc)
n=5;p=63;
x=(2*pi/p)*(0:p-1)';
y=mod(x+.44+.09*cos(3*x),2*pi);
PX=[ones(p,1),cos(x),sin(x),cos(2*x),sin(2*x)];
PY=[ones(p,1),cos(y),sin(y),cos(2*y),sin(2*y)];
PY(:,3)=PY(:,3)+.03*cos(7*x);   % reproducible contaminated snapshots
W=1+.4*cos(x).^2;
M=koopman_mpedmd(n,PX,PY,W,'metric');
U=koopman_mpedmd(n,PX,PY,W);
E=koopman_mpedmd(n,PX,PY,W,'edmd');
S=koopman_mpedmd(n,PX,PY,W,'spectral');
verifyTrue(tc,isreal(U));
verifyLessThan(tc,norm(U'*U-eye(n),'fro'),1e-10*n);
verifyLessThan(tc,norm(M.K'*M.G*M.K-M.G,'fro'),1e-9);
verifyLessThan(tc,norm(E-M.EDMD,'fro'),1e-12);
verifyGreaterThan(tc,norm(E'*M.G*E-M.G,'fro'),1e-4);
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
verifyEmpty(tc,S.finite.moments);
verifyLessThan(tc,abs(sum(M.weights)-1),1e-14);
end

function test_observable_weighted_spectral_measure_vs_trace(tc)
n=7;
M=koopman_mpedmd(n,[],[],[],'metric');
Q=M.Q;
v=complex(cos((1:n)'),sin(.31*(1:n)'));
[angles,weights]=um_koopman_cyclic_measure(Q,M.H,v);
verifyEqual(tc,numel(angles),n);
verifyEqual(tc,numel(weights),n);
verifyGreaterThanOrEqual(tc,min(weights),-1e-13);
verifyLessThan(tc,abs(sum(weights)-1),1e-12);
t=M.H*v;
for k=-3:3
    cyclic=sum(weights.*exp(1i*k*angles));
    exact=(t'*(Q^k)*t)/(t'*t);
    verifyLessThan(tc,abs(cyclic-exact),1e-10);
end
traceMoment=trace(Q)/n;
cyclicMoment=sum(weights.*exp(1i*angles));
verifyGreaterThan(tc,abs(traceMoment-cyclicMoment),1e-3);
end

function test_singular_cross_Gram_completion_and_bad_Gram(tc)
n=3;p=17;x=linspace(0,2*pi,p+1)';x(end)=[];
PX=[ones(p,1),cos(x),sin(x)];
W=ones(p,1)/p;
PY=zeros(p,n);
M=koopman_mpedmd(n,PX,PY,W,'metric');
verifyEqual(tc,M.cross_rank,0);
verifyTrue(tc,M.rank_deficient_completion);
verifyLessThan(tc,norm(M.Q'*M.Q-eye(n),'fro'),2e-12);
verifyLessThan(tc,norm(M.B,'fro'),2e-12);
verifyError(tc,@()koopman_mpedmd(n,PX,PX(:,1:2),W), ...
    'koopman_mpedmd:Dictionary');
bad=PX;bad(:,2)=PX(:,1);
verifyError(tc,@()koopman_mpedmd(n,bad,PX,W), ...
    'koopman_mpedmd:Gram');
verifyError(tc,@()koopman_mpedmd(n,PX,PX,zeros(p,1)), ...
    'koopman_mpedmd:Weights');
verifyError(tc,@()koopman_mpedmd(n,PX,PX,[ones(p-1,1);-1]), ...
    'koopman_mpedmd:Weights');
end

function test_custom_scaled_weights_and_real_data(tc)
n=3;p=50;x=linspace(0,2*pi,p+1)';x(end)=[];
PX=[ones(p,1),cos(x),sin(x)];
theta=.42;y=x+theta;
PY=[ones(p,1),cos(y),sin(y)];
A=koopman_mpedmd(n,PX,PY,ones(p,1),'metric');
B=koopman_mpedmd(n,PX,PY,200*ones(p,1),'metric');
verifyLessThan(tc,norm(A.Q-B.Q,'fro'),2e-12);
verifyLessThan(tc,norm(A.K-B.K,'fro'),2e-12);
verifyTrue(tc,isreal(A.Q));
verifyLessThan(tc,norm(A.Q'*A.Q-eye(n),'fro'),3e-12);
verifyLessThan(tc,norm(A.K'*A.G*A.K-A.G,'fro'),3e-12);
% Invariant real Fourier basis has rotation pair with phase theta.
verifyLessThan(tc,norm(sort(angle(eig(A.Q))) ...
   -[-theta;0;theta]),3e-10);
end

function test_input_modes_properties_and_rejection(tc)
[Q,p]=koopman_mpedmd();
verifyEmpty(tc,Q);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyTrue(tc,any(strcmp(p,'normal')));
verifyEqual(tc,size(koopman_mpedmd(4)),[4 4]);
verifyError(tc,@()koopman_mpedmd(0),'koopman_mpedmd:Order');
verifyError(tc,@()koopman_mpedmd(4,rand(10,4),[]), ...
    'koopman_mpedmd:Snapshots');
verifyError(tc,@()koopman_mpedmd(3,rand(4,3),rand(5,3)), ...
    'koopman_mpedmd:Dictionary');
verifyError(tc,@()koopman_mpedmd(4,[],[],[],'bad'), ...
    'koopman_mpedmd:Mode');
end
