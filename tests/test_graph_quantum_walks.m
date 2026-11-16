function tests=test_graph_quantum_walks
%TEST_GRAPH_QUANTUM_WALKS Independent spectra, local actions and NEPs.
tests=functiontests(localfunctions);
end

function test_szegedy_general_markov_spectrum(testCase)
cases={1;ones(3)/3; ...
    [.7 .3;.2 .8]; ...
    [0 1 0;0 0 1;1 0 0]; ...
    [.3 .4 .3;.2 .6 .2;.4 .1 .5]};
for j=1:numel(cases)
    P=cases{j};n=size(P,1);N=n*n;
    W=szegedy_walk(n,P);
    S=szegedy_walk(n,P,'spectral');
    expected=S.finite.eigenvalues();
    actual=eig(full(W));
    verifyEqual(testCase,numel(expected),N);
    verifyLessThan(testCase,norm(full(W'*W-eye(N)),'fro'),2e-11*N);
    verifyLessThan(testCase,hausdorff_spectrum(actual,expected),5e-10);
    for k=[0 1 2 3]
        verifyLessThan(testCase,abs(trace(W^k)/N-S.finite.moments(k)), ...
            2e-10);
    end
    O=szegedy_walk(n,P,'operator');
    check_operator(testCase,W,O);
end
end

function test_szegedy_large_matrix_free(testCase)
n=120;
O=szegedy_walk(n,[],'operator');
verifyEqual(testCase,O([],'size'),[n*n n*n]);
x=sin((1:n*n)')+1i*cos((1:n*n)'*.1);
y=O(x,'notransp');
verifyLessThan(testCase,norm(O(y,'ctransp')-x)/norm(x),5e-12);
S=szegedy_walk(n,[],'spectral');
verifyEqual(testCase,S.dimension,n*n);
verifyEqual(testCase,S.finite.kind,'exact');
end

function test_coined_cycle_closed_form(testCase)
for n=[3 4 7 10]
    A=sparse(1:n-1,2:n,1,n,n);A(1,n)=1;A=A+A';
    W=coined_graph_walk(n,A,'grover');
    verifyEqual(testCase,size(W),[2*n 2*n]);
    verifyLessThan(testCase,norm(full(W^n-speye(2*n)),'fro'),5e-11);
    z=exp(2i*pi*(0:n-1)'/n);
    verifyLessThan(testCase,hausdorff_spectrum(eig(full(W)),[z;z]),5e-10);
    check_operator(testCase,W,coined_graph_walk(n,A,'grover','operator'));
end
end

function test_coined_custom_coins(testCase)
A=[0 2 1 0;2 0 1 1;1 1 0 3;0 1 3 0];
[src,~,rev,d,~]=um_graph_arc_data(A);
coins=cell(1,4);
for j=1:4,coins{j}=eye(d(j));end
W=coined_graph_walk(4,A,coins);
S=sparse(rev,1:numel(src),1,numel(src),numel(src));
verifyEqual(testCase,W,S,'AbsTol',0);
verifyLessThan(testCase,norm(full(W^2-speye(numel(src))),'fro'),1e-13);
P=coined_graph_walk(4,A,'fourier','pencil');
verifyLessThan(testCase,norm(full(P.B\P.A-P.U),'fro'),2e-12);
verifyLessThan(testCase,norm(full(P.U'*P.U-speye(numel(src))),'fro'), ...
    2e-12);
check_operator(testCase,P.U,coined_graph_walk(4,A,'fourier','operator'));
end

function test_metric_graph_interval_secular_zeros(testCase)
ell=sqrt(2);
L=[0 ell;ell 0];
for k=[0 .4 pi/ell 2*pi/ell]
    W=metric_quantum_graph(2,L,k);
    desired=exp(1i*k*ell)*[0 1;1 0];
    verifyLessThan(testCase,norm(full(W)-desired,'fro'),1e-13);
    verifyLessThan(testCase,norm(full(W'*W-eye(2)),'fro'),1e-13);
    check_operator(testCase,W,metric_quantum_graph(2,L,k,[],'operator'));
end
F=metric_quantum_graph(2,L,.4,[],'secular');
verifyLessThan(testCase,abs(F.det(pi/ell)),5e-13);
verifyLessThan(testCase,abs(F.det(2*pi/ell)),5e-13);
verifyLessThan(testCase,norm(full(F.T(.4))- ...
    (eye(2)-full(metric_quantum_graph(2,L,.4))),'fro'),1e-13);
z=.3+.7i;
verifyLessThan(testCase,norm(full(F.T(z) - (eye(2)- ...
    exp(1i*z*ell)*[0 1;1 0])),'fro'),1e-12);
end

function test_metric_arbitrary_graph_local_formula(testCase)
L=[0 1.1 sqrt(2) 0; 1.1 0 .9 0; sqrt(2) .9 0 pi/3; 0 0 pi/3 0];
k=.73;
W=metric_quantum_graph(4,L,k);
[src,dst,~,d,~]=um_graph_arc_data(L);
N=numel(src);
ref=complex(zeros(N,N));
for from=1:N
    middle=dst(from);
    for to=find(src==middle)'
        backtrack=double(dst(to)==src(from));
        ref(to,from)=(2/d(middle)-backtrack)* ...
                    exp(1i*k*L(src(to),dst(to)));
    end
end
verifyLessThan(testCase,norm(full(W)-ref,'fro'),2e-12);
verifyLessThan(testCase,norm(full(W'*W-eye(N)),'fro'),2e-12);
check_operator(testCase,W,metric_quantum_graph(4,L,k,[],'operator'));
end

function test_modes_and_errors(testCase)
families={@szegedy_walk,@coined_graph_walk,@metric_quantum_graph};
for j=1:numel(families)
    [A,props]=families{j}();
    verifyEmpty(testCase,A);
    verifyTrue(testCase,ismember('unitary',props));
end
verifyError(testCase,@()szegedy_walk(3,ones(3)), ...
    'szegedy_walk:Stochastic');
verifyError(testCase,@()coined_graph_walk(3,eye(3)), ...
    'um_graph_arc_data:Adjacency');
verifyError(testCase,@()metric_quantum_graph(2,[0 1;1 0],.5, ...
    [],'pencil'),'metric_quantum_graph:Mode');
S=coined_graph_walk(4,[],[],'spectral');
verifyEqual(testCase,S.dimension,6);
T=metric_quantum_graph(3,[],.5,[],'spectral');
verifyEqual(testCase,T.dimension,6);
end

function check_operator(testCase,A,O)
N=size(A,1);
verifyEqual(testCase,O([],'size'),[N N]);
X=complex(sin((1:N)'*(1:2)),cos((1:N)'*(1:2)));
verifyLessThan(testCase,norm(O(X,'notransp')-A*X,'fro'),3e-11*N);
verifyLessThan(testCase,norm(O(X,'transp')-A'*X,'fro'),3e-11*N);
verifyLessThan(testCase,norm(O(O(X,'notransp'),'ctransp')-X,'fro'), ...
    3e-11*N);
end

function d=hausdorff_spectrum(a,b)
a=a(:);b=b(:);
d=max([max(arrayfun(@(z)min(abs(z-b)),a)), ...
       max(arrayfun(@(z)min(abs(z-a)),b))]);
end
