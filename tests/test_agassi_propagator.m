function tests=test_agassi_propagator
%TEST_AGASSI_PROPAGATOR Independent four-mode and number-sector checks.
tests=functiontests(localfunctions);
end

function test_hermiticity_unitarity_number_conservation(tc)
for j=[1 2]
    for pars=[1 .3 .4;1 0 .3;1 .3 0;0 .5 .2]'
        epsilon=pars(1);g=pars(2);V=pars(3);
        H=agassi_propagator(j,epsilon,g,V,.37,[],'hamiltonian');
        d=2^(4*j);
        verifyEqual(tc,size(H),[d d]);
        verifyTrue(tc,issparse(H));
        verifyLessThan(tc,norm(full(H-H'),'fro'),2e-12*d);
        s=uint32((0:d-1)');counts=zeros(d,1);
        for k=1:4*j,counts=counts+double(bitget(s,k));end
        Q=spdiags(counts,0,d,d);
        verifyLessThan(tc,norm(full(H*Q-Q*H),'fro'),3e-12*d);
        if j==1
            U=agassi_propagator(j,epsilon,g,V,.37);
            verifyLessThan(tc,norm(U'*U-eye(d),'fro'),3e-12*d);
        end
    end
end
end

function test_four_mode_published_pauli_spectrum(tc)
% Equations (14)-(17) in the four-qubit digital Agassi publication.
I=eye(2);z=diag([1 -1]);cp=[0 0;1 0];cm=cp';
epsilon=1.4;g=.3;V=.41;
H1=(epsilon-g)/4*(kron4(z,I,I,I)+kron4(I,z,I,I)) ...
   -(epsilon+g)/4*(kron4(I,I,z,I)+kron4(I,I,I,z));
H2=-g/4*(kron4(z,z,I,I)+kron4(I,I,z,z));
H3=-(g+V)*(kron4(cp,cp,cm,cm)+kron4(cm,cm,cp,cp));
paper=H1+H2+H3;
ours=full(agassi_propagator(1,epsilon,g,V,.4,[],'hamiltonian'));
% The occupation conventions differ by a unitary basis change and a
% constant -g/2 energy reference in this four-site specialization.
a=sort(real(eig(ours+g/2*eye(16))));
b=sort(real(eig(paper)));
verifyLessThan(tc,norm(a-b),2e-11);
end

function test_full_and_half_filling_sector(tc)
epsilon=1.2;g=.21;V=.37;t=.4;
H=full(agassi_propagator(1,epsilon,g,V,t,[],'hamiltonian'));
states=uint32((0:15)');counts=zeros(16,1);
for p=1:4,counts=counts+double(bitget(states,p));end
take=find(counts==2);
R=full(agassi_propagator(1,epsilon,g,V,t,2,'hamiltonian'));
verifyEqual(tc,size(R),[6 6]);
verifyLessThan(tc,norm(R-H(take,take),'fro'),5e-13);
U=agassi_propagator(1,epsilon,g,V,t,2);
verifyLessThan(tc,norm(U'*U-eye(6),'fro'),5e-13);
W=agassi_propagator(2,epsilon,g,V,t,4);
verifyEqual(tc,size(W),[70 70]);
verifyLessThan(tc,norm(W'*W-eye(70),'fro'),2e-11);
end

function test_free_spectral_phases(tc)
for j=[1 2]
  for particles={[],2*j}
    part=particles{1};
    W=agassi_propagator(j,1.3,0,0,.81,part);
    S=agassi_propagator(j,1.3,0,0,.81,part,'spectral');
    verifyEqual(tc,S.finite.kind,'exact');
    verifyEqual(tc,S.dimension,size(W,1));
    z=S.finite.eigenvalues();
    verifyLessThan(tc,norm(W-diag(z),'fro'),3e-12*size(W,1));
    verifyLessThan(tc,abs(trace(W^3)/size(W,1)-S.finite.moments(3)), ...
        2e-12);
  end
end
end

function test_metadata_and_validation(tc)
[e,p]=agassi_propagator();
verifyEmpty(tc,e);
verifyTrue(tc,any(strcmp(p,'unitary')));
S=agassi_propagator(1,[],[],[],[],[],'spectral');
verifyEqual(tc,S.family,'agassi_propagator');
verifyEqual(tc,S.finite.kind,'unknown');
T=agassi_propagator(1,[],[],[],0,[],'spectral');
verifyEqual(tc,T.finite.kind,'exact');
verifyError(tc,@()agassi_propagator(4),'MATLAB:notLessEqual');
verifyError(tc,@()agassi_propagator(1,[],[],[],[],5), ...
    'MATLAB:notLessEqual');
verifyError(tc,@()agassi_propagator(1,[],[],[],[],[], ...
    'operator'),'agassi_propagator:Mode');
end

function T=kron4(A,B,C,D)
T=kron(kron(kron(A,B),C),D);
end
