function tests=test_berry_wilson_loop
tests=functiontests(localfunctions);
end

function test_rank_one_latitude_and_quadratic_refinement(tc)
theta=1.1;
expected=exp(2i*pi*sin(theta/2)^2);
errors=zeros(1,3);
for j=1:3
    m=8*2^(j-1);
    frames=latitude_frames(m,theta);
    U=berry_wilson_loop(1,frames);
    D=berry_wilson_loop(1,frames,'data');
    S=berry_wilson_loop(1,frames,'spectral');
    verifyEqual(tc,size(U),[1 1]);
    verifyLessThan(tc,abs(abs(U)-1),2e-13);
    verifyEqual(tc,S.finite.kind,'exact');
    verifyLessThan(tc,abs(S.finite.eigenvalues()-U),2e-13);
    for k=-5:5
        verifyLessThan(tc,abs(S.finite.moments(k)-U^k),2e-13);
    end
    verifyLessThan(tc,abs(D.holonomy-U),2e-13);
    verifyLessThan(tc,abs(D.raw_product),1);
    errors(j)=abs(U-expected);
end
verifyGreaterThan(tc,errors(1)/errors(2),3.5);
verifyGreaterThan(tc,errors(2)/errors(3),3.5);
verifyLessThan(tc,errors(3),.005);
end

function test_nonabelian_gauge_covariance_and_raw_contractions(tc)
n=2;
data=berry_wilson_loop(n,[],'data');
verifyEqual(tc,data.ambient_dimension,4);
verifyEqual(tc,data.dimension,2);
verifyEqual(tc,data.samples,16);
verifyGreaterThan(tc,norm(data.links(:,:,1)*data.links(:,:,2)- ...
    data.links(:,:,2)*data.links(:,:,1),'fro'),.01);
verifyGreaterThan(tc,norm(data.raw_product'*data.raw_product-eye(n),'fro'),.1);
verifyLessThan(tc,norm(data.holonomy'*data.holonomy-eye(n),'fro'),3e-12);
% Reconstruct the same deterministic frames as the generator.
frames=default_rank_two_frames(16);
m=size(frames,3);
gauges=complex(zeros(n,n,m));
for j=1:m
    a=.12*j;b=.17*j;
    rot=[cos(a),-sin(a);sin(a),cos(a)];
    gauges(:,:,j)=diag(exp(1i*[b;-b/2]))*rot;
end
changed=frames;
for j=1:m,changed(:,:,j)=frames(:,:,j)*gauges(:,:,j);end
E=berry_wilson_loop(n,changed,'data');
verifyLessThan(tc,norm(E.holonomy- ...
    gauges(:,:,1)'*data.holonomy*gauges(:,:,1),'fro'),3e-11);
for j=1:m
    next=mod(j,m)+1;
    expected=gauges(:,:,j)'*data.links(:,:,j)*gauges(:,:,next);
    verifyLessThan(tc,norm(E.links(:,:,j)-expected,'fro'),2e-11);
end
for z=eig(data.holonomy).'
    verifyLessThan(tc,min(abs(eig(E.holonomy)-z)),2e-10);
end
S=berry_wilson_loop(n,changed,'spectral');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEmpty(tc,S.finite.eigenvalues);
end

function test_matrix_free_multi_rhs_and_default(tc)
for n=[1 2 3 5]
    U=berry_wilson_loop(n);
    F=berry_wilson_loop(n,[],'operator');
    X=randn(n,3)+1i*randn(n,3);
    verifyEqual(tc,size(U),[n n]);
    verifyEqual(tc,F([],'size'),[n n]);
    verifyLessThan(tc,norm(U'*U-eye(n),'fro'),6e-11*n);
    verifyLessThan(tc,norm(F(X,'notransp')-U*X,'fro'),8e-11*n);
    verifyLessThan(tc,norm(F(X,'transp')-U'*X,'fro'),8e-11*n);
    verifyLessThan(tc,norm(F(X,'ctransp')-U'*X,'fro'),8e-11*n);
    verifyLessThan(tc,norm(F(F(X,'notransp'),'ctransp')-X,'fro'),8e-11*n);
    links=berry_wilson_loop(n,[],'links');
    verifyEqual(tc,size(links,1),n);
    verifyEqual(tc,size(links,2),n);
    V=lattice_wilson_loop(n,links);
    verifyLessThan(tc,norm(V-U,'fro'),4e-12*n);
end
[a,p]=berry_wilson_loop();
verifyEmpty(tc,a);
verifyTrue(tc,any(strcmp(p,'unitary')));
end

function test_singular_overlaps_frame_validation_and_operator_errors(tc)
z=[1;0];w=[1e-11;sqrt(1-1e-22)];
frames=complex(zeros(2,1,2));
frames(:,:,1)=z;frames(:,:,2)=w;
verifyError(tc,@()berry_wilson_loop(1,frames), ...
    'berry_wilson_loop:SingularOverlap');
bad=frames;bad(:,:,2)=[.3;.5];
verifyError(tc,@()berry_wilson_loop(1,bad), ...
    'berry_wilson_loop:Orthonormal');
verifyError(tc,@()berry_wilson_loop(2,frames), ...
    'berry_wilson_loop:Frames');
verifyError(tc,@()berry_wilson_loop(2,[], 'bogus'), ...
    'berry_wilson_loop:Mode');
F=berry_wilson_loop(2,[],'operator');
verifyError(tc,@()F(ones(3,1),'notransp'), ...
    'um_wilson_apply:Dimension');
verifyError(tc,@()F(ones(2,1),'bogus'), ...
    'um_wilson_apply:Flag');
S=berry_wilson_loop(1,[],'spectral');
verifyError(tc,@()S.finite.moments(1.5), ...
    'berry_wilson_loop:Moments');
end

function frames=latitude_frames(m,theta)
frames=complex(zeros(2,1,m));
for j=1:m
    t=2*pi*(j-1)/m;
    frames(:,1,j)=[cos(theta/2);exp(1i*t)*sin(theta/2)];
end
end

function frames=default_rank_two_frames(m)
d=4;n=2;
Q0=fft(eye(d))/sqrt(d);
Q0=Q0(:,1:n);
B=complex(zeros(d));
for j=1:d-1
    B(j,j+1)=.38+.07i*j;
    B(j+1,j)=-conj(B(j,j+1));
end
frames=complex(zeros(d,n,m));
for j=1:m
    t=2*pi*(j-1)/m;
    frames(:,:,j)=diag(exp(1i*(0:d-1)*t))*expm(sin(t)*B)*Q0;
end
end
