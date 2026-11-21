function tests=test_flag_connection_expmv
%TEST_FLAG_CONNECTION_EXPMV Independent dense and geometric flag checks.
tests=functiontests(localfunctions);
end

function test_exact_dense_connection_match_complex_multirhs(tc)
cases={[1 1 1],[2 2],[1 2 2],[2 2 2],[1 2 3],[1 1 1 1],[3 2 2]};
for j=1:numel(cases)
    blocks=cases{j};n=sum(blocks);
    d=(n*n-sum(blocks.^2))/2;
    ref=real_flag_parallel_transport(blocks,[],.73,'model');
    for t=[-.7 0 .91 5.3]
        Q=real_flag_parallel_transport(blocks, ...
            ref.horizontal_geodesic_generator,t,'matrix');
        O=real_flag_parallel_transport(blocks, ...
            ref.horizontal_geodesic_generator,t,'operator');
        X=reshape(sin(.033*(1:d*3)),d,3) ...
           +1i*reshape(cos(.023*(1:d*3)),d,3);
        verifyEqual(tc,O([],'size'),[d d]);
        verifyLessThan(tc,norm(O(X,'notransp')-Q*X,'fro'), ...
            4e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(O(X,'ctransp')-Q'*X,'fro'), ...
            4e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(O(X,'transp')-Q.'*X,'fro'), ...
            4e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X, ...
            'fro'),5e-11*max(1,norm(X,'fro')));
        verifyLessThan(tc,abs(norm(O(X,'notransp'),'fro') ...
            -norm(X,'fro')),4e-11*max(1,norm(X,'fro')));
    end
end
end

function test_independent_lie_projector_tangent_generator(tc)
for blocks={[1 1 2],[2 1 2],[1 1 1 1]}
    sizes=blocks{1};
    M=real_flag_parallel_transport(sizes,[],.2,'model');
    n=sum(sizes);d=M.dimension;
    K=M.horizontal_geodesic_generator;
    O=real_flag_parallel_transport(sizes,K,1e-5,'operator');
    % Difference quotient is an independent projection of
    % [K,F_ij] constructed in full Lie algebra.
    for q=1:min(d,6)
        E=zeros(n);
        E(M.basis_row(q),M.basis_col(q))=1;
        E(M.basis_col(q),M.basis_row(q))=-1;
        C=K*E-E*K;
        expected=-.5*C(sub2ind([n n],M.basis_row,M.basis_col));
        x=zeros(d,1);x(q)=1;
        numerical=(O(x,'notransp')-x)/1e-5;
        verifyLessThan(tc,norm(numerical-expected), ...
            2e-5*max(1,norm(K,'fro')^2));
        verifyLessThan(tc,norm(expected-M.connection_generator(:,q)), ...
            3e-12*d);
    end
end
end

function test_large_flag_without_d_square_storage(tc)
% A two-block Grassmann moving-frame connection is identically zero.
sizes=[128 128];d=16384;
O=real_flag_parallel_transport(sizes,[],.7,'operator');
X=reshape(cos(.05*(1:2*d)),d,2);
verifyEqual(tc,O([],'size'),[d d]);
verifyEqual(tc,O(X,'notransp'),X);
verifyEqual(tc,O(X,'ctransp'),X);
verifyError(tc,@()real_flag_parallel_transport(sizes,[],.7), ...
    'real_flag_parallel_transport:DenseOrder');
% Nontrivial three-block example beyond dense n and d limits.
sizes=[2 3 91];n=sum(sizes);d=(n*n-sum(sizes.^2))/2;
K=sparse(n,n);
K(1,3)=.17;K(3,1)=-.17;
K(3,6)=-.29;K(6,3)=.29;
K(2,7)=.13;K(7,2)=-.13;
O=real_flag_parallel_transport(sizes,K,.84,'operator');
X=reshape(sin(.022*(1:d*2)),d,2) ...
    +1i*reshape(cos(.037*(1:d*2)),d,2);
Y=O(X,'notransp');
verifyLessThan(tc,norm(O(Y,'ctransp')-X,'fro'), ...
    1e-10*max(1,norm(X,'fro')));
verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')), ...
    1e-10*max(1,norm(X,'fro')));
verifyGreaterThan(tc,norm(Y-X,'fro'),.05);
end

function test_invalid_operator_arguments(tc)
O=real_flag_parallel_transport([2 3 1],[],.4,'operator');
verifyError(tc,@()O(zeros(10,1),'notransp'), ...
    'um_flag_connection_expmv:Dimension');
verifyError(tc,@()O(zeros(11,1),'wrong'), ...
    'um_flag_connection_expmv:Flag');
verifyError(tc,@()real_flag_parallel_transport([257 1],[],.5, ...
    'operator'),'real_flag_parallel_transport:Order');
end
