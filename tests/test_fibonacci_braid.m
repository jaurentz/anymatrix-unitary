function tests=test_fibonacci_braid
%TEST_FIBONACCI_BRAID Fusion paths, braid identities and spectral data.
tests=functiontests(localfunctions);
end

function test_fibonacci_sector_dimensions_and_orthogonality(tc)
fib=ones(1,20);
for j=3:20,fib(j)=fib(j-1)+fib(j-2);end
for n=3:9
    for charge={'tau','vacuum'}
        label=charge{1};
        data=fibonacci_braid(n,label,[],'generators');
        if strcmp(label,'tau'),expected=fib(n);else,expected=fib(n-1);end
        verifyEqual(tc,data.dimension,expected);
        verifyEqual(tc,size(data.fusion_paths),[expected n]);
        verifyEqual(tc,data.fusion_paths(:,1),ones(expected,1));
        for k=1:numel(data.braid_generators)
            G=data.braid_generators{k};
            verifyTrue(tc,issparse(G));
            verifyLessThanOrEqual(tc,nnz(G),2*expected);
            verifyLessThan(tc,norm(full(G'*G-speye(expected)),'fro'), ...
                5e-12*max(1,expected));
        end
    end
end
end

function test_yang_baxter_and_distant_braids(tc)
for n=3:9
    for charge={'tau','vacuum'}
        for chirality=[1 -1]
            data=fibonacci_braid(n,charge{1},[],'generators',chirality);
            g=data.braid_generators;
            for k=1:n-2
                a=g{k};b=g{k+1};
                verifyLessThan(tc,norm(full(a*b*a-b*a*b),'fro'), ...
                    3e-11*data.dimension);
            end
            for k=1:n-1
                for j=k+2:n-1
                    verifyLessThan(tc,norm(full(g{k}*g{j}-g{j}*g{k}),'fro'), ...
                        3e-11*data.dimension);
                end
            end
        end
    end
end
end

function test_three_anyon_F_R_and_four_anyon_vacuum(tc)
phi=(1+sqrt(5))/2;
F=[1/phi,1/sqrt(phi);1/sqrt(phi),-1/phi];
r=[exp(-4i*pi/5);exp(3i*pi/5)];
D=fibonacci_braid(3,'tau',[],'generators');
verifyEqual(tc,D.fusion_paths,[1 0 1;1 1 1]);
verifyLessThan(tc,norm(D.F-F,'fro'),1e-13);
verifyLessThan(tc,norm(D.braid_generators{1}-diag(r),'fro'),1e-13);
verifyLessThan(tc,norm(D.braid_generators{2}-F*diag(r)*F,'fro'),1e-13);
D4=fibonacci_braid(4,'vacuum',[],'generators');
verifyEqual(tc,D4.dimension,2);
verifyLessThan(tc,norm(D4.braid_generators{1}-D4.braid_generators{3},'fro'),1e-13);
verifyLessThan(tc,norm(D4.braid_generators{2}-F*diag(r)*F,'fro'),1e-13);
end

function test_exact_braid_generator_phase_multiplicities(tc)
for n=3:8
    for charge={'tau','vacuum'}
        for chirality=[1 -1]
            d=fibonacci_braid(n,charge{1},[],'generators',chirality);
            for k=1:n-1
                for sign=[1 -1]
                    S=fibonacci_braid(n,charge{1},sign*k,'spectral',chirality);
                    lam=eig(full(d.braid_generators{k}));
                    if sign<0,lam=conj(lam);end
                    ref=S.finite.eigenvalues();
                    verifyEqual(tc,S.finite.kind,'exact');
                    verifyEqual(tc,numel(ref),d.dimension);
                    verifyEqual(tc,sum(S.finite.atoms.weights),1,'AbsTol',1e-14);
                    for t=-2:4
                        verifyLessThan(tc,abs(mean(lam.^t)-S.finite.moments(t)), ...
                            2e-11);
                    end
                    for z=unique(ref(:)).'
                        verifyEqual(tc,sum(abs(ref-z)<1e-9), ...
                            sum(abs(lam-z)<1e-9));
                    end
                end
            end
        end
    end
end
end

function test_braid_word_operator_inverse_and_chirality(tc)
n=8;charge='tau';word=[1 2 -3 4 6 -7 2];
U=fibonacci_braid(n,charge,word);
O=fibonacci_braid(n,charge,word,'operator');
d=size(U,1);
X=randn(d,3)+1i*randn(d,3);
verifyEqual(tc,O([],'size'),[d d]);
verifyLessThan(tc,norm(U'*U-eye(d),'fro'),1e-10);
verifyLessThan(tc,norm(O(X,'notransp')-U*X,'fro'),1e-10);
verifyLessThan(tc,norm(O(X,'ctransp')-U'*X,'fro'),1e-10);
verifyLessThan(tc,norm(O(X,'transp')-U'*X,'fro'),1e-10);
verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X,'fro'),1e-10);
inverse=fibonacci_braid(n,charge,-fliplr(word));
verifyLessThan(tc,norm(inverse-U','fro'),1e-10);
other=fibonacci_braid(n,charge,word,'matrix',-1);
verifyLessThan(tc,norm(other-conj(U),'fro'),1e-10);
unknown=fibonacci_braid(n,charge,word,'spectral');
verifyEqual(tc,unknown.finite.kind,'unknown');
identity=fibonacci_braid(n,charge,[],'spectral');
verifyEqual(tc,identity.finite.kind,'exact');
verifyLessThan(tc,abs(identity.finite.moments(4)-1),1e-14);
end

function test_large_sparse_and_validation(tc)
d=fibonacci_braid(16,'tau',[],'generators');
verifyEqual(tc,d.dimension,987);
word=[1 2 -3 8 -10 14];
O=fibonacci_braid(16,'tau',word,'operator');
X=complex(sin((1:987)'*.15),cos((1:987)'*.31));
verifyEqual(tc,O([],'size'),[987 987]);
verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X), ...
    1e-10*norm(X));
verifyError(tc,@()fibonacci_braid(16,'tau',word), ...
    'fibonacci_braid:ExplicitOrder');
[A,props]=fibonacci_braid();
verifyEmpty(tc,A);
verifyTrue(tc,any(strcmp(props,'unitary')));
verifyError(tc,@()fibonacci_braid(3,'notacharge'), ...
    'fibonacci_braid:Charge');
verifyError(tc,@()fibonacci_braid(3,'tau',0), ...
    'fibonacci_braid:Word');
verifyError(tc,@()fibonacci_braid(3,'tau',3), ...
    'fibonacci_braid:Word');
verifyError(tc,@()fibonacci_braid(3,'tau',1,'bad'), ...
    'fibonacci_braid:Mode');
verifyError(tc,@()fibonacci_braid(3,'tau',1,'matrix',0), ...
    'fibonacci_braid:Chirality');
verifyError(tc,@()O(ones(986,1),'notransp'), ...
    'um_fibonacci_braid_apply:Dimension');
verifyError(tc,@()O(X,'invalid'), ...
    'um_fibonacci_braid_apply:Flag');
end
