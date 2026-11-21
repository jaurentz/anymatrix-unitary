function tests=test_gauss_large_order_stability
%TEST_GAUSS_LARGE_ORDER_STABILITY Sturm signs, log Christoffel weights.
tests=functiontests(localfunctions);
end

function test_hermite_log_weights_and_function_scaling(tc)
for n=[1 2 6 24 64 192]
    Q=gauss_hermite_transform(n);
    S=gauss_hermite_transform(n,'quadrature');
    verifyEqual(tc,S.family,'gauss_hermite_transform');
    verifyEqual(tc,numel(S.nodes),n);
    verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),3e-10*n);
    verifyLessThan(tc,S.orthogonality_residual,3e-10*n);
    verifyLessThan(tc,S.diagonalization_residual,3e-10*n);
    verifyLessThan(tc,abs(sum(S.normalized_weights)-1),2e-8);
    verifyTrue(tc,all(isfinite(S.log_normalized_weights)));
    verifyTrue(tc,all(diff(S.nodes)>0) || n==1);
    verifyGreaterThanOrEqual(tc,min(Q(:,1)),-5e-13);
    if n>1
        verifyLessThan(tc,max(abs(S.nodes+flipud(S.nodes))),1e-9*n);
    end
    reliable=find(S.normalized_weights>1e-10);
    verifyLessThan(tc,max(abs(log(Q(reliable,1).^2) - ...
        S.log_normalized_weights(reliable))),3e-7);
    % Independently evaluated ground Hermite function
    phi0=pi^(-1/4)*exp(-.5*S.nodes.^2);
    if n<=64
        verifyLessThan(tc,norm(Q(:,1)./exp(S.log_function_scaling) ...
            -phi0),2e-10*n);
    end
    Sphase=gauss_hermite_transform(n,'spectral');
    verifyEqual(tc,Sphase.finite.kind,'unknown');
end
end

function test_laguerre_log_weights_near_minus_one(tc)
for alpha=[-.999 -.9 0 .5 3 16]
    for n=[1 5 24 64 160]
        Q=gauss_laguerre_transform(n,alpha);
        S=gauss_laguerre_transform(n,alpha,'quadrature');
        verifyEqual(tc,S.family,'gauss_laguerre_transform');
        verifyEqual(tc,S.alpha,alpha);
        verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),4e-10*n);
        verifyLessThan(tc,S.orthogonality_residual,4e-10*n);
        verifyLessThan(tc,S.diagonalization_residual,2e-8*n);
        verifyGreaterThan(tc,min(S.nodes),0);
        verifyLessThan(tc,abs(sum(S.normalized_weights)-1),3e-7);
        verifyTrue(tc,all(isfinite(S.log_normalized_weights)));
        verifyTrue(tc,all(isfinite(S.log_function_scaling)));
        verifyGreaterThanOrEqual(tc,min(Q(:,1)),-3e-12);
        reliable=find(S.normalized_weights>1e-9);
        verifyLessThan(tc,max(abs(log(Q(reliable,1).^2)- ...
            S.log_normalized_weights(reliable))),2e-6);
        if n<=64
            phi0=exp(.5*alpha*log(S.nodes)-.5*S.nodes ...
                -.5*gammaln(alpha+1));
            verifyLessThan(tc,norm(Q(:,1)./ ...
                exp(S.log_function_scaling)-phi0),2e-8*n);
        end
    end
    Sphase=gauss_laguerre_transform(7,alpha,'spectral');
    verifyEqual(tc,Sphase.finite.kind,'unknown');
end
end

function test_extreme_tails_and_orthogonal_actions(tc)
n=384;
for which=1:2
    if which==1
        Q=gauss_hermite_transform(n);
        S=gauss_hermite_transform(n,'quadrature');
    else
        Q=gauss_laguerre_transform(n,.37);
        S=gauss_laguerre_transform(n,.37,'quadrature');
    end
    verifyTrue(tc,all(isfinite(S.nodes)));
    verifyTrue(tc,all(isfinite(S.log_normalized_weights)));
    verifyLessThan(tc,min(S.log_normalized_weights),-700);
    verifyLessThan(tc,abs(sum(S.normalized_weights)-1),8e-7);
    verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),1e-8);
    verifyGreaterThanOrEqual(tc,min(Q(:,1)),-8e-13);
    X=randn(n,3)+1i*randn(n,3);
    verifyLessThan(tc,norm(Q'*(Q*X)-X,'fro'),2e-8);
    verifyLessThan(tc,norm(Q*(Q'*X)-X,'fro'),2e-8);
    verifyTrue(tc,all(isfinite(S.log_function_scaling)));
end
end

function test_invalid_modes_and_gauss_sturm_sign(tc)
verifyError(tc,@()gauss_hermite_transform(7,'operator'), ...
    'gauss_hermite_transform:Mode');
verifyError(tc,@()gauss_laguerre_transform(7,1,'operator'), ...
    'gauss_laguerre_transform:Mode');
verifyError(tc,@()gauss_laguerre_transform(7,-1), ...
    'MATLAB:notGreater');
for n=[2 11 57]
    d=linspace(-3,4,n)';
    e=linspace(.2,3,n-1)';
    Q=um_gauss_basis(d,e);
    verifyLessThan(tc,norm(Q'*Q-eye(n),'fro'),8e-11*n);
    verifyGreaterThanOrEqual(tc,min(Q(:,1)),-5e-14);
    last=Q(:,end);parity=(-1).^(n-(1:n))';
    reliable=abs(last)>1e-10;
    verifyTrue(tc,all(last(reliable).*parity(reliable)>0));
end
end
