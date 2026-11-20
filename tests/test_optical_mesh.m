function tests=test_optical_mesh
tests=functiontests(localfunctions);
end

function test_reck_clements_independent_core_products(tc)
for n=[1 2 3 4 5 6]
    m=n*(n-1)/2;
    theta=.15+.9*(1:m)'/max(1,m+1);
    phi=(1:m)'*.37;
    ext=(0:n-1)'*.21;
    for layout={'reck','clements'}
        name=layout{1};
        A=optical_mesh(n,name,theta,phi,ext);
        verifyEqual(tc,size(A),[n n]);
        verifyLessThan(tc,norm(A'*A-eye(n),'fro'),3e-12*n);
        verifyLessThan(tc,abs(abs(det(A))-1),3e-12*n);
        verifyLessThan(tc,abs(det(A)-exp(1i*(sum(phi)+sum(ext)))), ...
            3e-11*n);
        p=independent_schedule(n,name);
        H=diag(exp(1i*ext));
        Q=eye(n);
        for t=1:m
            j=p(t);c=cos(theta(t));s=sin(theta(t));
            G=eye(n);
            G(j:j+1,j:j+1)=[exp(1i*phi(t))*c,-s; ...
                exp(1i*phi(t))*s,c];
            Q=G*Q;
        end
        verifyLessThan(tc,norm(A-H*Q,'fro'),3e-12*n);
        F=optical_mesh(n,name,theta,phi,ext,'operator');
        X=randn(n,3)+1i*randn(n,3);
        verifyEqual(tc,F([],'size'),[n n]);
        verifyLessThan(tc,norm(F(X,'notransp')-A*X,'fro'),3e-12*n);
        verifyLessThan(tc,norm(F(X,'transp')-A'*X,'fro'),3e-12*n);
        verifyLessThan(tc,norm(F(X,'ctransp')-A'*X,'fro'),3e-12*n);
        verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'), ...
            3e-12*n);
        S=optical_mesh(n,name,theta,phi,ext,'spectral');
        if n==1,verifyEqual(tc,S.finite.kind,'exact');
        else,verifyEqual(tc,S.finite.kind,'unknown');end
    end
end
end

function test_balanced_coupler_and_routing(tc)
M=optical_mesh(2,'reck',pi/4,0,[0;0]);
verifyEqual(tc,M,[1 -1;1 1]/sqrt(2),'AbsTol',2e-14);
for n=[2 3 4 5]
    m=n*(n-1)/2;
    U=optical_mesh(n,'clements',(pi/2)*ones(m,1),zeros(m,1), ...
        zeros(n,1));
    verifyLessThan(tc,norm(U'*U-eye(n),'fro'),3e-13*n);
    verifyLessThan(tc,norm(abs(U)*ones(n,1)-ones(n,1)),3e-13*n);
    verifyLessThan(tc,norm(abs(U)'*ones(n,1)-ones(n,1)),3e-13*n);
    verifyLessThan(tc,norm(abs(U).^2-abs(U),'fro'),3e-13*n);
end
end

function test_phase_spectral_moments_and_sensitivity(tc)
n=5;m=n*(n-1)/2;
theta=zeros(m,1);phi=.3*(1:m)';ext=.2*(0:n-1)';
S=optical_mesh(n,'reck',theta,phi,ext,'spectral');
U=optical_mesh(n,'reck',theta,phi,ext);
verifyEqual(tc,S.finite.kind,'exact');
verifyEqual(tc,numel(S.finite.eigenvalues()),n);
for k=-4:4
    verifyLessThan(tc,abs(S.finite.moments(k)-trace(U^k)/n),4e-12);
end
theta(:)=.23;
step=1e-6;
v=randn(n,1)+1i*randn(n,1);
a=optical_mesh(n,'reck',theta,phi,ext);
tp=theta;tp(3)=tp(3)+step;
tm=theta;tm(3)=tm(3)-step;
fd=(optical_mesh(n,'reck',tp,phi,ext)- ...
    optical_mesh(n,'reck',tm,phi,ext))/(2*step);
verifyLessThan(tc,norm(fd*v,'fro'),2);
verifyGreaterThan(tc,norm(fd*v,'fro'),1e-5);
verifyLessThan(tc,norm(a'*a-eye(n),'fro'),3e-12);
end

function test_invalid_inputs(tc)
verifyError(tc,@()optical_mesh(3,'bad'),'optical_mesh:Layout');
verifyError(tc,@()optical_mesh(4,'reck',zeros(2,1)), ...
    'optical_mesh:Theta');
verifyError(tc,@()optical_mesh(3,'reck',[.1;.2;pi]), ...
    'optical_mesh:Theta');
verifyError(tc,@()optical_mesh(3,'reck',zeros(3,1),zeros(2,1)), ...
    'optical_mesh:Phi');
verifyError(tc,@()optical_mesh(3,'reck',zeros(3,1),zeros(3,1),[0 1]), ...
    'optical_mesh:External');
verifyError(tc,@()optical_mesh(3,'reck',[],[],[],'wrong'), ...
    'optical_mesh:Mode');
F=optical_mesh(4,'clements',[],[],[],'operator');
verifyError(tc,@()F(ones(3,2),'notransp'), ...
    'um_optical_mesh_apply:Dimension');
verifyError(tc,@()F(ones(4,2),'bad'), ...
    'um_optical_mesh_apply:Flag');
end

function pairs=independent_schedule(n,layout)
pairs=[];
switch layout
    case 'reck'
        for col=1:n-1
            pairs=[pairs, n-1:-1:col]; %#ok<AGROW>
        end
    case 'clements'
        for layer=1:n
            pairs=[pairs, (1+mod(layer-1,2)):2:n-1]; %#ok<AGROW>
        end
end
verify_count=numel(pairs);
if verify_count~=n*(n-1)/2
    error('test_optical_mesh:Internal','Unexpected number of couplers.');
end
end
