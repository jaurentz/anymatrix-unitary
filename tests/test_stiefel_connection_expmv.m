function tests=test_stiefel_connection_expmv
%TEST_STIEFEL_CONNECTION_EXPMV Matrix-free canonical tangent transport.
tests=functiontests(localfunctions);
end

function test_dense_reference_and_complex_forward_adjoint(tc)
cases=[3 1;4 2;7 3;10 4;15 2;48 4;7 7];
for h=1:size(cases,1)
    n=cases(h,1);p=cases(h,2);d=n*p-p*(p+1)/2;
    A=zeros(p);
    for i=1:p
        for j=i+1:p
            A(i,j)=.35*sin(.43*i+.27*j);
            A(j,i)=-A(i,j);
        end
    end
    B=reshape(.24*sin((1:(n-p)*p)*.37),n-p,p);
    t=(-1)^h * (1.2+.13*h);
    Q=stiefel_parallel_transport(n,p,A,B,t);
    O=stiefel_parallel_transport(n,p,A,B,t,'operator');
    X=reshape(sin(.014*(1:3*d)),d,3) ...
        +1i*reshape(cos(.023*(1:3*d)),d,3);
    verifyEqual(tc,O([],'size'),[d d]);
    verifyLessThan(tc,norm(O(X,'notransp')-Q*X,'fro'), ...
        6e-11*max(1,norm(X,'fro')));
    verifyLessThan(tc,norm(O(X,'ctransp')-Q'*X,'fro'), ...
        6e-11*max(1,norm(X,'fro')));
    verifyLessThan(tc,norm(O(X,'transp')-Q.'*X,'fro'), ...
        6e-11*max(1,norm(X,'fro')));
    verifyLessThan(tc,norm(O(O(X,'notransp'),'ctransp')-X, ...
        'fro'),5e-11*max(1,norm(X,'fro')));
    verifyLessThan(tc,abs(norm(O(X,'notransp'),'fro') ...
        -norm(X,'fro')),4e-11*max(1,norm(X,'fro')));
end
end

function test_step_scaling_and_short_long_time_behavior(tc)
n=9;p=3;d=n*p-p*(p+1)/2;
A=[0 -.5 .13;.5 0 -.27;-.13 .27 0];
B=reshape(.4*cos(.23*(1:(n-p)*p)),n-p,p);
X=reshape(cos(.11*(1:2*d)),d,2);
for t=[0 1e-5 .13 2.6 11.4 -8.7]
    D=stiefel_parallel_transport(n,p,A,B,t);
    O=stiefel_parallel_transport(n,p,A,B,t,'operator');
    verifyLessThan(tc,norm(O(X,'notransp')-D*X,'fro'), ...
        1e-10*max(1,norm(X,'fro')));
end
end

function test_large_n_small_p_subquadratic_storage_and_isometry(tc)
% Exceed the dense n<=64,d<=256 caps. The operator creates no d-by-d
% connection or n-by-n SO(n) generator.
for dims=[128 2;256 3;512 2;1024 2]'
    n=dims(1);p=dims(2);d=n*p-p*(p+1)/2;
    A=zeros(p);
    if p>1,A(1,2)=.27;A(2,1)=-.27;end
    B=reshape(.032*sin(.01*(1:(n-p)*p)),n-p,p);
    O=stiefel_parallel_transport(n,p,A,B,.73,'operator');
    X=reshape(sin(.07*(1:d*2)),d,2);
    Y=O(X,'notransp');
    verifyEqual(tc,size(Y),size(X));
    verifyEqual(tc,O([],'size'),[d d]);
    verifyLessThan(tc,norm(O(Y,'ctransp')-X,'fro'), ...
        2e-10*max(1,norm(X,'fro')));
    verifyLessThan(tc,abs(norm(Y,'fro')-norm(X,'fro')), ...
        2e-10*max(1,norm(X,'fro')));
    if n>64
        verifyError(tc,@()stiefel_parallel_transport(n,p,A,B,.73), ...
            'stiefel_parallel_transport:Order');
    end
end
end

function test_sphere_group_limit_and_invalid_flags(tc)
n=150;p=1;B=.7*cos(.1*(1:n-1)');
O=stiefel_parallel_transport(n,p,0,B,.7,'operator');
X=reshape(cos(.2*(1:2*(n-1))),n-1,2);
verifyEqual(tc,O(X,'notransp'),X);
verifyEqual(tc,O(X,'ctransp'),X);
O=stiefel_parallel_transport(9,2,[],[],0,'operator');
X=reshape(1:30,15,2);
verifyEqual(tc,O(X,'notransp'),X);
O=stiefel_parallel_transport(9,2,[],[],.5,'operator');
verifyError(tc,@()O(zeros(14,1),'notransp'), ...
    'um_stiefel_connection_expmv:Dimension');
verifyError(tc,@()O(zeros(15,1),'bad'), ...
    'um_stiefel_connection_expmv:Flag');
verifyError(tc,@()stiefel_parallel_transport(110,40,[],[],.4,'operator'), ...
    'stiefel_parallel_transport:OperatorSize');
end
