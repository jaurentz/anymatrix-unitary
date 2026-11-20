function tests=test_resonant_scattering
%TEST_RESONANT_SCATTERING Flux conservation, resonances and delays.
tests=functiontests(localfunctions);
end

function test_on_shell_unitarity_and_reciprocity(tc)
for n=[1 2 3 7 16]
    for E=[-3 -.4 0 .27 2.1]
        U=resonant_scattering(n,E);
        verifyEqual(tc,size(U),[n n]);
        verifyLessThan(tc,norm(U'*U-eye(n),'fro'),2e-11*n);
        verifyLessThan(tc,norm(U-U.','fro'),2e-11*n);
        X=reshape(sin((1:3*n)*.39),n,3) ...
              +1i*reshape(cos((1:3*n)*.26),n,3);
        verifyLessThan(tc,norm(U*X,'fro')-norm(X,'fro'), ...
            2e-11*max(1,norm(X,'fro')));
    end
end
U=resonant_scattering(3,0.1);
verifyGreaterThan(tc,norm(U-diag(diag(U)),'fro'),.01);
end

function test_wigner_smith_derivative_and_positive_delay(tc)
for n=[1 2 5 9]
    for E=[-1.1 .13 1.7]
        U=resonant_scattering(n,E);
        D=resonant_scattering(n,E,[],[],'derivative');
        Q=resonant_scattering(n,E,[],[],'time_delay');
        h=2e-6;
        F=(resonant_scattering(n,E+h)-resonant_scattering(n,E-h))/(2*h);
        verifyLessThan(tc,norm(F-D,'fro'),2e-7*max(1,norm(D,'fro')));
        verifyLessThan(tc,norm(Q+1i*U'*D,'fro'),2e-10*n);
        verifyLessThan(tc,norm(Q-Q','fro'),1e-12*n);
        verifyGreaterThan(tc,min(eig((Q+Q')/2)),0);
    end
end
H=[-.4 .13i;-.13i .6];
W=[.65+.1i, .13, .05;-.09i,.7,.18i];
E=.37;
U=resonant_scattering(3,E,H,W);
Q=resonant_scattering(3,E,H,W,'time_delay');
D=resonant_scattering(3,E,H,W,'derivative');
verifyLessThan(tc,norm(U'*U-eye(3),'fro'),4e-12);
verifyLessThan(tc,norm(Q+1i*U'*D,'fro'),4e-11);
verifyGreaterThan(tc,norm(U-U.','fro'),1e-3);
verifyGreaterThanOrEqual(tc,min(eig(Q)),-1e-12);
end

function test_resonance_poles_complex_continuation_and_internal_evp(tc)
model=resonant_scattering(4,.1,[],[],'model');
verifyEqual(tc,model.dimension,4);
verifyEqual(tc,model.internal_order,4);
verifyLessThan(tc,norm(model.H-model.H','fro'),1e-13);
verifyLessThan(tc,norm(model.Gamma-model.W*model.W','fro'),1e-13);
verifyTrue(tc,all(imag(model.poles)<-1e-6));
verifyLessThan(tc,norm(model.S(.17)-resonant_scattering(4,.17),'fro'),2e-12);
verifyLessThan(tc,norm(model.time_delay(.17) ...
    -resonant_scattering(4,.17,[],[],'time_delay'),'fro'),2e-12);
verifyLessThan(tc,norm(sort(eig(model.effective_hamiltonian)) ...
    -sort(model.poles)),2e-12);
H=.8; W=sqrt(.3);
d=resonant_scattering(1,.5,H,W,'model');
verifyLessThan(tc,abs(d.poles-(.8-.3i)),1e-14);
verifyLessThan(tc,abs(eig(H)-.8),1e-14);
near=d.S(d.poles+1e-4);
verifyGreaterThan(tc,abs(near),1000);
verifyLessThan(tc,abs(d.S(.8)) - 1,1e-12);
verifyError(tc,@()d.time_delay(.8+.1i), ...
    'um_resonant_scattering_eval:Energy');
end

function test_decoupled_exact_phases_and_lorentzian_delay(tc)
epsi=[-.7;.4;1.2];gamma=[.2;.4;.6];
H=diag(epsi);W=diag(sqrt(gamma));
for E=[-1 0 .17 .7 1.8]
    S=resonant_scattering(3,E,H,W,'spectral');
    A=resonant_scattering(3,E,H,W);
    expected=(E-epsi-1i*gamma)./(E-epsi+1i*gamma);
    verifyEqual(tc,S.finite.kind,'exact');
    verifyEqual(tc,S.dimension,3);
    verifyLessThan(tc,norm(A-diag(expected),'fro'),4e-12);
    verifyLessThan(tc,norm(S.finite.eigenvalues()-expected),2e-12);
    verifyEqual(tc,sum(S.finite.atoms.weights),1,'AbsTol',1e-13);
    for k=-4:5
        verifyLessThan(tc,abs(S.finite.moments(k)-mean(expected.^k)), ...
            4e-12);
    end
    Q=resonant_scattering(3,E,H,W,'time_delay');
    expectedQ=2*gamma./((E-epsi).^2+gamma.^2);
    verifyLessThan(tc,norm(Q-diag(expectedQ),'fro'),3e-12);
end
% Unwrap the single-channel eigenphase through its resonance.
grid=linspace(-1.2,1.2,901);
s=zeros(size(grid));delay=s;
for j=1:numel(grid)
    s(j)=resonant_scattering(1,grid(j),0,sqrt(.3));
    delay(j)=resonant_scattering(1,grid(j),0,sqrt(.3),'time_delay');
end
phase=unwrap(angle(s));
slope=(phase(3:end)-phase(1:end-2))./(grid(3:end)-grid(1:end-2));
verifyLessThan(tc,max(abs(slope-delay(2:end-1))),2.5e-4);
verifyTrue(tc,all(diff(phase)>-1e-9));
end

function test_rectangular_channels_and_contractive_observed_block(tc)
H=diag([-.3 .6]);
W=[sqrt(.3) 0 0 0;0 sqrt(.25) sqrt(.25) 0];
E=.2;n=4;
U=resonant_scattering(n,E,H,W);
Q=resonant_scattering(n,E,H,W,'time_delay');
S=resonant_scattering(n,E,H,W,'spectral');
verifyEqual(tc,S.finite.kind,'exact');
verifyEqual(tc,numel(S.finite.eigenvalues()),n);
verifyLessThan(tc,norm(U'*U-eye(n),'fro'),3e-12);
verifyLessThan(tc,norm(U-U.','fro'),3e-12);
verifyLessThan(tc,norm(U*U'-eye(n),'fro'),3e-12);
verifyLessThan(tc,abs(sum(abs(eig(U)-1)<1e-10)-2),1e-12);
verifyLessThan(tc,norm(Q-Q','fro'),1e-13);
verifyEqual(tc,sum(eig(Q)>1e-9),2);
% Omitting outgoing/incoming channel 3 loses flux into that channel.
keep=[1 2 4];
B=U(keep,keep);
loss=eye(3)-B'*B;
verifyGreaterThanOrEqual(tc,min(eig((loss+loss')/2)),-1e-12);
verifyGreaterThan(tc,norm(loss,'fro'),1e-3);
verifyLessThanOrEqual(tc,max(eig(B'*B)),1+1e-12);
end

function test_public_api_and_invalid_models(tc)
[Z,p]=resonant_scattering();
verifyEmpty(tc,Z);
verifyTrue(tc,any(strcmp(p,'unitary')));
verifyTrue(tc,any(strcmp(p,'square')));
verifyEqual(tc,resonant_scattering(1,0,0,1),-1,'AbsTol',1e-14);
S=resonant_scattering(1,.1,[],[],'spectral');
verifyEqual(tc,S.finite.kind,'exact');
T=resonant_scattering(3,.1,[],[],'spectral');
verifyEqual(tc,T.finite.kind,'unknown');
verifyEmpty(tc,T.finite.eigenvalues);
verifyError(tc,@()resonant_scattering(0), ...
    'resonant_scattering:Channels');
verifyError(tc,@()resonant_scattering(2,0,[1 1;0 1],eye(2)), ...
    'resonant_scattering:Hamiltonian');
verifyError(tc,@()resonant_scattering(2,0,eye(2), ...
    [1 0;1 0]),'resonant_scattering:Coupling');
verifyError(tc,@()resonant_scattering(1,0,eye(2), ...
    [1;1]),'resonant_scattering:Channels');
verifyError(tc,@()resonant_scattering(2,0,eye(2),1), ...
    'resonant_scattering:CouplingShape');
verifyError(tc,@()resonant_scattering(2,1i), ...
    'resonant_scattering:Energy');
verifyError(tc,@()resonant_scattering(2,0,[],[],'bogus'), ...
    'resonant_scattering:Mode');
end
