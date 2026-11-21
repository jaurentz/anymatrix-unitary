function Y=um_flat_torus_spin_dirac_apply(X,flag,nx,ny,t,mass,twist,px,py,E)
%UM_FLAT_TORUS_SPIN_DIRAC_APPLY Exact finite modal 2D spinor FFT action.
%   The FFT diagonalizes the finite twisted shift lattice exactly in
%   mathematical arithmetic. All output modes share the same action.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_flat_torus_spin_dirac_apply:Flag','FLAG must be scalar text.');
end
flag=lower(char(flag));
n=nx*ny;
if strcmp(flag,'size')
    Y=[2*n 2*n];return
end
if ~any(strcmp(flag,{'notransp','ctransp','transp'}))
    error('um_flat_torus_spin_dirac_apply:Flag','Invalid action flag.');
end
if ~(isnumeric(X)&&ismatrix(X)&&size(X,1)==2*n)
    error('um_flat_torus_spin_dirac_apply:Dimension', ...
        'X must have 2*NX*NY rows.');
end
if strcmp(flag,'transp')
    % U.' X = conjugate(U' conjugate(X)); distinct from U' in general.
    Y=conj(um_flat_torus_spin_dirac_apply(conj(X),'ctransp', ...
        nx,ny,t,mass,twist,px,py,E));
    return
end
if strcmp(flag,'ctransp'),t=-t;end
c=cos(t*E);
s=zeros(size(E));
small=(E<1e-14);
s(small)=t;
s(~small)=sin(t*E(~small))./E(~small);
u11=c-1i*mass*s;
u22=c+1i*mass*s;
u12=-1i*s.*(px-1i*py);
u21=-1i*s.*(px+1i*py);
[x,y]=ndgrid((0:nx-1)',(0:ny-1)');
g=exp(1i*(twist(1)*x/nx+twist(2)*y/ny));
Y=complex(zeros(size(X)));
for rhs=1:size(X,2)
    a=reshape(X(1:2:end,rhs),nx,ny);
    b=reshape(X(2:2:end,rhs),nx,ny);
    fa=fft2(conj(g).*a);
    fb=fft2(conj(g).*b);
    ya=g.*ifft2(u11.*fa+u12.*fb);
    yb=g.*ifft2(u21.*fa+u22.*fb);
    Y(1:2:end,rhs)=ya(:);
    Y(2:2:end,rhs)=yb(:);
end
end
