function tests=test_real_flag_parallel_transport
%TEST_REAL_FLAG_PARALLEL_TRANSPORT Intrinsic quotient LC validation.
tests=functiontests(localfunctions);
end

function test_independent_projector_metric_connection_and_isometry(tc)
% Compare Lie-bracket generator with the independently derived
% induced metric connection on the complete tuple of flag projectors.
cases={[1 1 1],[2 2],[1 2 2],[1 1 2 1],[2 2 2]};
for h=1:numel(cases)
    blocks=cases{h};
    n=sum(blocks);
    d=(n^2-sum(blocks.^2))/2;
    M=real_flag_parallel_transport(blocks,[],.71,'model');
    Q=M.U;S=M.connection_generator;K=M.horizontal_geodesic_generator;
    verifyEqual(tc,size(Q),[d d]);
    verifyLessThan(tc,norm(Q'*Q-eye(d),'fro'),4e-11*d);
    verifyLessThan(tc,norm(S+S','fro'),3e-12*d);
    verifyLessThan(tc,norm(M.moving_frame'*M.moving_frame ...
        -eye(n),'fro'),3e-12*n);
    F=cell(d,1);
    for j=1:d
        A=zeros(n);i=M.basis_row(j);k=M.basis_col(j);
        A(i,k)=1;A(k,i)=-1;F{j}=A;
    end
    Gram=zeros(d);connection=zeros(d);
    for i=1:d
        for j=1:d
            for a=1:numel(blocks)
                P=M.initial_projectors{a};
                Ei=F{i}*P-P*F{i};
                Ej=F{j}*P-P*F{j};
                Ejdot=K*Ej-Ej*K;
                Gram(i,j)=Gram(i,j)+trace(Ei'*Ej)/4;
                connection(i,j)=connection(i,j)+trace(Ei'*Ejdot)/4;
            end
        end
    end
    verifyLessThan(tc,norm(Gram-eye(d),'fro'),3e-12*d);
    verifyLessThan(tc,norm(connection+S,'fro'),3e-12*d);
    verifyLessThan(tc,norm(Q-expm(-.71*connection),'fro'), ...
        3e-11*d);
    for a=1:numel(blocks)
        P=M.initial_projectors{a};Pt=M.endpoint_projectors{a};
        verifyLessThan(tc,norm(Pt-Pt','fro'),3e-12*n);
        verifyLessThan(tc,norm(Pt*Pt-Pt,'fro'),3e-12*n);
        verifyLessThan(tc,norm(Pt-M.moving_frame*P ...
            *M.moving_frame','fro'),3e-12*n);
        verifyLessThan(tc,abs(trace(Pt)-blocks(a)),3e-12*n);
    end
    % Initial geodesic horizontal velocity stays parallel.
    v=zeros(d,1);
    for i=1:d
        v(i)=K(M.basis_row(i),M.basis_col(i));
    end
    verifyLessThan(tc,norm(Q*v-v),3e-11*max(1,norm(v)));
    Bwd=real_flag_parallel_transport(blocks,K,-.71);
    verifyLessThan(tc,norm(Bwd*Q-eye(d),'fro'),3e-11*d);
end
end

function test_grassmann_two_block_identity_and_nontrivial_three_block(tc)
for blocks={[1 2],[2 3],[3 4]}
    sizes=blocks{1};
    n=sum(sizes);d=prod(sizes);
    U=real_flag_parallel_transport(sizes,[],.9);
    verifyLessThan(tc,norm(U-eye(d),'fro'),2e-12*d);
    S=real_flag_parallel_transport(sizes,[],.9,'spectral');
    verifyEqual(tc,S.dimension,d);
    verifyEqual(tc,S.finite.kind,'exact');
    verifyEqual(tc,S.finite.eigenvalues(),ones(d,1));
    verifyEqual(tc,S.finite.moments([0 1 2 -3]),ones(1,4));
end
A=real_flag_parallel_transport([1 2 2],[],.83);
verifyGreaterThan(tc,norm(A-eye(size(A)),'fro'),.1);
S=real_flag_parallel_transport([1 2 2],[],.83,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
end

function test_full_flag_biinvariant_group_half_adjoint(tc)
for n=[3 4 5]
    blocks=ones(1,n);t=.61;
    M=real_flag_parallel_transport(blocks,[],t,'model');
    K=M.horizontal_geodesic_generator;
    H=expm(-t*K/2);
    d=n*(n-1)/2;
    B=zeros(d);
    for j=1:d
        F=zeros(n);
        a=M.basis_row(j);b=M.basis_col(j);
        F(a,b)=1;F(b,a)=-1;
        transported=H*F*H';
        for i=1:d
            B(i,j)=transported(M.basis_row(i),M.basis_col(i));
        end
    end
    verifyLessThan(tc,norm(B-M.U,'fro'),3e-12*d);
end
end

function test_invalid_and_reproducible(tc)
[U,props]=real_flag_parallel_transport();
verifyEmpty(tc,U);
verifyTrue(tc,any(strcmp(props,'orthogonal')));
A=real_flag_parallel_transport([1 2 1],[],1);
B=real_flag_parallel_transport([1 2 1],[],1);
verifyEqual(tc,A,B,'AbsTol',0);
verifyError(tc,@()real_flag_parallel_transport([1 0 2]), ...
    'real_flag_parallel_transport:Blocks');
verifyError(tc,@()real_flag_parallel_transport([3]), ...
    'real_flag_parallel_transport:Blocks');
verifyError(tc,@()real_flag_parallel_transport([1 1]), ...
    'real_flag_parallel_transport:Order');
verifyError(tc,@()real_flag_parallel_transport([5 30]), ...
    'real_flag_parallel_transport:Order');
n=4;
verifyError(tc,@()real_flag_parallel_transport([2 2],eye(n)), ...
    'real_flag_parallel_transport:Skew');
K=zeros(n);K(1,2)=.3;K(2,1)=-.3;
verifyError(tc,@()real_flag_parallel_transport([2 2],K), ...
    'real_flag_parallel_transport:Horizontal');
verifyError(tc,@()real_flag_parallel_transport([2 2],[],NaN), ...
    'real_flag_parallel_transport:Time');
verifyError(tc,@()real_flag_parallel_transport([2 2],[],1, ...
    'operator'),'real_flag_parallel_transport:Mode');
end
