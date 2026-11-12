function tests = test_core_representations
%TEST_CORE_REPRESENTATIONS Unitary, structural and spectral invariants.
tests=functiontests(localfunctions);
end

function test_all_shapes_and_parameters(testCase)
param=[0 0; -0.5 -0.5; 0.5 0.5; 0.3 1.6; -0.6 0.4; 2.4 -0.7];
forms={'hessenberg','inverse_hessenberg','cmv','twisted'};
for ia=1:size(param,1)
    a=param(ia,1); b=param(ia,2);
    mu=szego_jacobi_measure(a,b);
    for n=[1 2 3 5 8 12]
        for beta=[1,-1,exp(1i*0.37)]
            H=szego_matrix(mu,n,'hessenberg',beta);
            for ir=1:numel(forms)
                U=szego_matrix(mu,n,forms{ir},beta);
                verifyEqual(testCase,size(U),[n n]);
                verifyTrue(testCase,issparse(U));
                verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),2e-12*max(n,1));
                if beta==1 || beta==-1
                    verifyTrue(testCase,isreal(U));
                end
                % Power traces independently certify equality of multisets.
                % For these small n, testing moments up to n suffices.
                X=speye(n); Y=speye(n);
                for k=1:n
                    X=X*H; Y=Y*U;
                    verifyLessThan(testCase,abs(trace(X)-trace(Y)),1e-8*n);
                end
                [z,w]=szego_quadrature(U);
                verifyLessThan(testCase,abs(sum(w)-1),1e-11);
                verifyGreaterThanOrEqual(testCase,min(w),0);
                reference=js_jacobi_cosine_moments(a,b,min(n-1,5));
                for k=1:min(n-1,5)
                    verifyLessThan(testCase,abs(sum(w.*z.^k)-reference(k+1)),1e-8);
                end
                if strcmp(forms{ir},'cmv')
                    verifyLessThanOrEqual(testCase,nnz(U),5*n);
                elseif strcmp(forms{ir},'hessenberg')
                    verifyLessThan(testCase,norm(tril(full(U),-2),'fro'),1e-12*n);
                elseif strcmp(forms{ir},'inverse_hessenberg')
                    verifyLessThan(testCase,norm(triu(full(U),2),'fro'),1e-12*n);
                end
            end
        end
    end
end
end

function test_arbitrary_snakes(testCase)
n=9; alpha=szego_verblunsky(szego_jacobi_measure(.25,.75),n-1);
twists=logical([1 1 0 1 0 0 1]);
[U,order]=js_core_product(alpha,exp(1i*.2),'twisted',twists);
H=js_core_product(alpha,exp(1i*.2),'hessenberg');
verifyEqual(testCase,sort(order),1:n-1);
verifyLessThan(testCase,norm(full(U'*U-speye(n)),'fro'),1e-11);
for k=1:n
    verifyLessThan(testCase,abs(trace(U^k)-trace(H^k)),1e-8*n);
    e=sparse(1,1,1,n,1);
    verifyLessThan(testCase,abs(e'*U^k*e-e'*H^k*e),1e-8);
end
end

function test_terminal_phase_for_order_one(testCase)
for beta=[1,-1,exp(1i*.31)]
    for form={'cmv','hessenberg','inverse_hessenberg','twisted'}
        U=js_core_product(zeros(0,1),beta,form{1});
        verifyEqual(testCase,full(U),conj(beta),'AbsTol',1e-14);
    end
end
end

function test_expected_errors(testCase)
verifyError(testCase,@()js_core_product(1,1,'cmv'),'js_core_product:Alpha');
verifyError(testCase,@()js_core_product(0,0.5,'cmv'),'js_core_product:Beta');
verifyError(testCase,@()js_core_product(0,1,'invalid'),'js_core_product:Representation');
verifyError(testCase,@()js_core_product([.1;.2;.3],1,'twisted',[1]),'js_core_product:Twists');
end
