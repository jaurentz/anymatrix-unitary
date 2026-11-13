function tests=test_dirac_operator
%TEST_DIRAC_OPERATOR Independent checks of Dirac matrix-free interface.
tests=functiontests(localfunctions);
end

function test_forward_adjoint_and_spectrum(tc)
s=rng;cleanup=onCleanup(@()rng(s)); %#ok<NASGU>
rng(317,'twister');
for n=[1 2 3 12]
    for mass=[0 .35 .8 1]
        A=dirac_automaton(n,mass);
        O=dirac_automaton(n,mass,'operator');
        verifyEqual(tc,O([],'size'),[2*n 2*n]);
        X=randn(2*n,3)+1i*randn(2*n,3);
        verifyLessThan(tc,norm(O(X,'notransp')-A*X,'fro'),1e-11);
        verifyLessThan(tc,norm(O(X,'transp')-A'*X,'fro'),1e-11);
        verifyLessThan(tc,norm(O(X,'ctransp')-A'*X,'fro'),1e-11);
        verifyLessThan(tc,norm(O(O(X,'notransp'),'transp')-X,'fro'),1e-11);
        omega=acos(sqrt(1-mass^2)*cos(2*pi*(0:n-1)'/n));
        exact=sort(angle(exp(1i*[omega;-omega])));
        computed=sort(angle(eig(full(A))));
        verifyLessThan(tc,max(abs(exact-computed)),1e-10);
    end
end
end

function test_bad_input(tc)
O=dirac_automaton(5,.4,'operator');
verifyError(tc,@()O(zeros(9,1),'notransp'), ...
    'um_dirac_apply:Dimension');
verifyError(tc,@()O(zeros(10,1),'nonsense'), ...
    'um_dirac_apply:Flag');
verifyError(tc,@()dirac_automaton(5,.4,'unknown'), ...
    'um_is_operator:Mode');
end
