function tests=test_wigner_skew_expmv
%TEST_WIGNER_SKEW_EXPMV Arbitrary-angle no-dense Wigner actions (#41).
tests=functiontests(localfunctions);
end

function test_generic_half_integer_and_integer_dense_euler_reference(tc)
for ell=[.5 1 2 3.5 5 10 16]
    n=round(2*ell)+1;
    X=reshape(sin(.39*(1:3*n)),n,3) ...
        +1i*reshape(cos(.23*(1:3*n)),n,3);
    for beta=[.31 -.79 1.9 pi+.037]
        a=.27;g=-.43;
        D=wigner_rotation(ell,a,beta,g);
        F=wigner_rotation(ell,a,beta,g,'operator');
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(F(X,'notransp')-D*X,'fro'), ...
            7e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(F(X,'transp')-D.'*X,'fro'), ...
            7e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(F(X,'ctransp')-D'*X,'fro'), ...
            7e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X, ...
            'fro'),1e-10*max(1,norm(X,'fro')));
        verifyLessThan(tc,abs(norm(F(X,'notransp'),'fro') ...
            -norm(X,'fro')),7e-11*norm(X,'fro'));
    end
end
end

function test_large_spin_no_dense_exponential(tc)
% A full 1026-by-1026 Wigner D is not assembled anywhere in this test.
ell=512.5;n=round(2*ell)+1;
j=(1:n)';
X=[sin(.031*j)+1i*cos(.044*j),cos(.017*j)+1i*sin(.025*j)];
F=wigner_rotation(ell,.29,.083,-.36,'operator');
verifyEqual(tc,F([],'size'),[n n]);
Y=F(X,'notransp');
verifyLessThan(tc,norm(F(Y,'ctransp')-X,'fro'), ...
    1e-9*norm(X,'fro'));
verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')), ...
    1e-9*norm(X,'fro'));
% Independent group law about y: D_y(b1)*D_y(b2)=D_y(b1+b2).
F1=wigner_rotation(ell,0,.043,0,'operator');
F2=wigner_rotation(ell,0,.052,0,'operator');
F12=wigner_rotation(ell,0,.095,0,'operator');
verifyLessThan(tc,norm(F1(F2(X,'notransp'),'notransp') ...
    -F12(X,'notransp'),'fro'),2e-9*norm(X,'fro'));
end

function test_complex_direct_sum_transpose_semantics(tc)
for L=[1 2 5 9]
    n=L^2;
    D=spherical_harmonic_rotation(L,.47,.61,-.24);
    F=spherical_harmonic_rotation(L,.47,.61,-.24,'operator');
    X=reshape(sin(.13*(1:3*n)),n,3) ...
      +1i*reshape(cos(.07*(1:3*n)),n,3);
    verifyEqual(tc,F([],'size'),[n n]);
    verifyLessThan(tc,norm(F(X,'notransp')-D*X,'fro'), ...
        8e-11*max(1,norm(X,'fro')));
    verifyLessThan(tc,norm(F(X,'transp')-D.'*X,'fro'), ...
        8e-11*max(1,norm(X,'fro')));
    verifyLessThan(tc,norm(F(X,'ctransp')-D'*X,'fro'), ...
        8e-11*max(1,norm(X,'fro')));
    if L>1
        verifyGreaterThan(tc,norm(F(X,'transp') ...
            -F(X,'ctransp'),'fro'),.01);
    end
end
end

function test_real_block_matrix_free_and_sparse_basis(tc)
for L=[1 2 5 10]
    n=L^2;
    D=real_spherical_harmonic_rotation(L,.21,.67,-.43);
    F=real_spherical_harmonic_rotation(L,.21,.67,-.43,'operator');
    X=reshape(sin(.17*(1:2*n)),n,2);
    Z=complex(X,reshape(cos(.11*(1:2*n)),n,2));
    verifyEqual(tc,F([],'size'),[n n]);
    verifyTrue(tc,isreal(F(X,'notransp')));
    verifyLessThan(tc,norm(F(X,'notransp')-D*X,'fro'), ...
        1e-10*max(1,norm(X,'fro')));
    verifyLessThan(tc,norm(F(Z,'notransp')-D*Z,'fro'), ...
        1e-10*max(1,norm(Z,'fro')));
    verifyLessThan(tc,norm(F(Z,'ctransp')-D'*Z,'fro'), ...
        1e-10*max(1,norm(Z,'fro')));
    verifyLessThan(tc,norm(F(Z,'transp')-D.'*Z,'fro'), ...
        1e-10*max(1,norm(Z,'fro')));
    for ell=0:L-1
        C=um_wigner_real_basis(ell);
        verifyTrue(tc,issparse(C));
        verifyEqual(tc,nnz(C),4*ell+1);
        verifyLessThan(tc,norm(C'*C-speye(2*ell+1),'fro'), ...
            2e-12*(2*ell+1));
    end
end
end

function test_high_bandlimit_no_cached_dense_blocks(tc)
% L=48 means 2304-dimensional full real/complex SO(3) rotations.
% Run an adjoint roundtrip without any O(L^4) dense construction.
L=48;n=L^2;
X=sin((1:n)'*.014)+1i*cos((1:n)'*.017);
for kind=1:2
    if kind==1
        F=spherical_harmonic_rotation(L,.21,.17,-.44,'operator');
    else
        F=real_spherical_harmonic_rotation(L,.21,.17,-.44,'operator');
    end
    Y=F(X,'notransp');
    verifyEqual(tc,F([],'size'),[n n]);
    verifyLessThan(tc,norm(F(Y,'ctransp')-X), ...
        5e-9*norm(X));
end
end

function test_invalid_operator_inputs(tc)
F=wigner_rotation(3,.2,.47,-.1,'operator');
verifyError(tc,@()F(zeros(6,1),'notransp'), ...
    'um_wigner_skew_expmv:Dimension');
verifyError(tc,@()F(ones(7,1),'wrong'), ...
    'um_wigner_skew_expmv:Flag');
end
