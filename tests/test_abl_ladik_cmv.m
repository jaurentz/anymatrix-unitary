function tests=test_abl_ladik_cmv
%TEST_ABL_LADIK_CMV Independent finite Lax and exact weight flow invariants.
tests=functiontests(localfunctions);
end

function test_exact_unitary_spectrum_and_weight_flow(tc)
for n=[4 6 8 12]
    z=exp(2i*pi*((0:n-1)'+.5)/n);
    w0=(1:n)'.^1.3;w0=w0/sum(w0);
    for k=[1 2 3]
        for t=[-.25 0 .4]
            A=full(abl_ladik_cmv(n,t,z,w0,k));
            S=abl_ladik_cmv(n,t,z,w0,k,'spectral');
            expected=w0.*exp(-2*t*real(z.^k));
            expected=expected/sum(expected);
            verifyEqual(tc,S.family,'abl_ladik_cmv');
            verifyEqual(tc,S.finite.kind,'exact');
            verifyLessThan(tc,norm(A'*A-eye(n),'fro'),5e-10*n);
            verifyLessThan(tc,norm(A-diag(diag(A)),'fro'),n); %#ok<NASGU>
            verifyLessThan(tc,max(max(abs(triu(A,3)))),2e-12*n);
            lam=eig(A);
            for j=1:n
                verifyLessThan(tc,min(abs(lam-z(j))),2e-8);
            end
            verifyLessThan(tc,norm(S.cyclic.weights-expected),3e-13);
            for p=-3:3
                verifyLessThan(tc,abs(trace(A^abs(p))/n- ...
                    sum(z.^abs(p))/n),8e-9*n);
                verifyLessThan(tc,abs(S.finite.moments(p)- ...
                    sum(z.^p)/n),1e-12);
                verifyLessThan(tc,abs(S.cyclic.moments(p)- ...
                    sum(expected.*z.^p)),1e-12);
            end
            par=abl_ladik_cmv(n,t,z,w0,k,'parameters');
            verifyLessThan(tc,abs(par.beta+1),2e-12);
            verifyLessThan(tc,abs(sum(par.weight_derivative)),3e-13);
            verifyLessThan(tc,max(abs(par.alpha)),1);
        end
    end
end
end

function test_lax_equation_against_independent_finite_differences(tc)
% The integrable Lax commutator is independently computed from C(t).
for n=[4 6 10]
    z=exp(2*pi*1i*((0:n-1)'+.5)/n);
    w=(1:n)'.^1.1;w=w/sum(w);
    for harmonic=[1 2 3]
        t=.17;h=2e-6;
        F=abl_ladik_cmv(n,t,z,w,harmonic,'lax');
        Cp=full(abl_ladik_cmv(n,t+h,z,w,harmonic));
        Cm=full(abl_ladik_cmv(n,t-h,z,w,harmonic));
        deriv=(Cp-Cm)/(2*h);
        verifyLessThan(tc,norm(F.B+F.B','fro'),1e-12*n);
        verifyLessThan(tc,norm(deriv-F.commutator,'fro'),4e-7*n);
        verifyLessThan(tc,norm(F.commutator-F.matrix*F.B+ ...
            F.B*F.matrix,'fro'),3e-13*n);
        verifyLessThan(tc,abs(sum(F.weight_derivative)),1e-12);
        verifyLessThan(tc,abs(trace(F.commutator)),5e-12*n);
    end
end
end

function test_real_schur_dmKdV_reduction(tc)
n=8;z=exp(2*pi*1i*((0:n-1)'+.5)/n);
w=(1:n)';w=w+flipud(w);w=w/sum(w);
t=.15;h=1e-6;
A=full(abl_ladik_cmv(n,t,z,w));
P=abl_ladik_cmv(n,t,z,w,1,'parameters');
verifyLessThan(tc,norm(imag(A),'fro'),1e-11*n);
verifyLessThan(tc,max(abs(imag(P.alpha))),1e-11);
Pa=abl_ladik_cmv(n,t+h,z,w,1,'parameters');
Pb=abl_ladik_cmv(n,t-h,z,w,1,'parameters');
adot=(Pa.alpha-Pb.alpha)/(2*h);
aa=[-1;P.alpha;-1];
% The real Schur flow is alpha_dot=(1-alpha^2)*(alpha_(j-1)-alpha_(j+1)).
rhs=(1-abs(P.alpha).^2).*(aa(1:end-2)-aa(3:end));
verifyLessThan(tc,norm(adot-rhs),1e-7*n);
end

function test_operator_and_initial_condition(tc)
n=12;z=exp(2*pi*1i*((0:n-1)'+.5)/n);
w=(1:n)';w=w/sum(w);
A=full(abl_ladik_cmv(n,.21,z,w,2));
F=abl_ladik_cmv(n,.21,z,w,2,'operator');
X=randn(n,3)+1i*randn(n,3);
verifyEqual(tc,F([],'size'),[n n]);
verifyLessThan(tc,norm(F(X,'notransp')-A*X,'fro'),3e-10*n);
verifyLessThan(tc,norm(F(X,'ctransp')-A'*X,'fro'),3e-10*n);
verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),4e-10*n);
A0=abl_ladik_cmv(n,0,z,w);
B0=inverse_opuc(n,z,w);
verifyLessThan(tc,norm(full(A0-B0),'fro'),2e-11*n);
end

function test_invalid_inputs(tc)
verifyError(tc,@()abl_ladik_cmv(5), 'abl_ladik_cmv:Order');
n=4;z=exp(2i*pi*((0:n-1)'+.5)/n);
verifyError(tc,@()abl_ladik_cmv(n,0,z*exp(.1i)), ...
    'abl_ladik_cmv:Boundary');
verifyError(tc,@()abl_ladik_cmv(n,0,z,[0;1;2;3]), ...
    'inverse_opuc:Weights');
verifyError(tc,@()abl_ladik_cmv(n,0,z,[1;1;1;1],0), ...
    'MATLAB:expectedPositive');
verifyError(tc,@()abl_ladik_cmv(n,0,z,[1;1;1;1],1,'bad'), ...
    'abl_ladik_cmv:Mode');
end
