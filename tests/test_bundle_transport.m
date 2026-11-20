function tests=test_bundle_transport
tests=functiontests(localfunctions);
end
function test_connection_cycle_and_degenerate_cycles(tc)
for n=[1 2 5 8]
    angles=.2+.1*cos(2*pi*(1:n)'/(n+1));
    R=zeros(2,2,n);
    for j=1:n
        a=angles(j);R(:,:,j)=[cos(a) -sin(a);sin(a) cos(a)];
    end
    p=[2:n 1];
    U=bundle_transport(n,p,R);
    verifyTrue(tc,issparse(U));
    verifyLessThan(tc,norm(full(U-connection_cycle(n,angles)),'fro'),1e-12);
    verifyLessThan(tc,norm(full(U'*U-eye(2*n)),'fro'),2e-12);
    compare(tc,U,bundle_transport(n,p,R,'spectral'));
end
p=[2 1 4 3];R=repmat(eye(2),1,1,4);
S=bundle_transport(4,p,R,'spectral');
verifyEqual(tc,numel(S.finite.eigenvalues()),8);
verifyLessThan(tc,abs(S.finite.moments(1)),1e-12);
verifyLessThan(tc,abs(S.finite.moments(2)-1),1e-12);
end
function test_nonabelian_gauge_and_operators(tc)
n=6;d=3;p=[2 3 1 5 4 6]';R=zeros(d,d,n);
for j=1:n
    A=zeros(d);A(1,2)=.11*j;A(2,1)=-.11*j;
    A(2,3)=.07*j;A(3,2)=-.07*j;
    R(:,:,j)=expm(A);
end
R(:,:,6)=diag(exp(1i*[.3 -.2 .8]));
U=bundle_transport(n,p,R);
verifyLessThan(tc,norm(full(U'*U-eye(n*d)),'fro'),2e-12);
S=bundle_transport(n,p,R,'spectral');compare(tc,U,S);
X=randn(n*d,3)+1i*randn(n*d,3);
F=bundle_transport(n,p,R,'operator');
verifyEqual(tc,F([],'size'),[n*d n*d]);
verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),2e-12);
verifyLessThan(tc,norm(F(X,'transp')-U'*X,'fro'),2e-12);
verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),2e-12);
verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),2e-12);
G=zeros(d,d,n);D=eye(n*d);
for j=1:n
    G(:,:,j)=diag(exp(1i*.1*j*(1:d)));
    idx=(j-1)*d+(1:d);D(idx,idx)=G(:,:,j);
end
Rg=zeros(d,d,n);
for j=1:n,Rg(:,:,j)=G(:,:,p(j))*R(:,:,j)*G(:,:,j)';end
Ug=bundle_transport(n,p,Rg);
verifyLessThan(tc,norm(full(Ug-D*U*D'),'fro'),2e-12);
compare(tc,Ug,bundle_transport(n,p,Rg,'spectral'));
end
function test_input_errors(tc)
L=repmat(eye(2),1,1,4);
verifyError(tc,@()bundle_transport(4,[1 1 3 4],L), ...
    'bundle_transport:Permutation');
L(:,:,2)=2*eye(2);
verifyError(tc,@()bundle_transport(4,[2 3 4 1],L), ...
    'bundle_transport:Unitary');
verifyError(tc,@()bundle_transport(4,[],ones(2,3,4)), ...
    'bundle_transport:Links');
F=bundle_transport(4,[],[],'operator');
verifyError(tc,@()F(ones(7,1),'notransp'), ...
    'um_bundle_transport_apply:Dimension');
verifyError(tc,@()F(ones(8,1),'bad'), ...
    'um_bundle_transport_apply:Flag');
end
function compare(tc,U,S)
verifyEqual(tc,S.finite.kind,'exact');
z=S.finite.eigenvalues();e=eig(full(U));
verifyEqual(tc,numel(z),numel(e));
for j=1:numel(e)
    verifyLessThan(tc,min(abs(z-e(j))),2e-10);
end
for k=[-5 -2 -1 0 1 2 3 6]
    verifyLessThan(tc,abs(S.finite.moments(k)-trace(full(U)^k)/size(U,1)), ...
        3e-10);
end
end
