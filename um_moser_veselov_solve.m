function F=um_moser_veselov_solve(h,J,M)
%UM_MOSER_VESELOV_SOLVE Near-identity special-orthogonal implicit solver.
%   Solve F'*J-J*F=h*M using X=F'=(I-K)\(I+K), K=-K'.
%   Newton residual R=X*J-J*X'-h*M in strict upper coordinates.
%   Analytic differential dX=(I-K)\(dK*(I+X)).
%   K initial guess solves (K*J+J*K)=h*M/2.
%   If Newton is unable to find the locally admissible branch, ERROR;
%   there is no universal solvability or uniqueness claim.
n=size(J,1);
I=eye(n);
target=h*M;
[V,D]=eig(J);
d=real(diag(D));
T=V'*target*V;
Khat=T./(2*(d+d.'));
K=V*Khat*V';
K=(K-K')/2;
[ii,jj]=find(triu(ones(n),1));
nd=numel(ii);
idx=sub2ind([n n],ii,jj);
tol=4e-12*max([1,norm(target,'fro'),norm(J,'fro')]);
solved=false;
for iter=1:35
    X=(I-K)\(I+K);
    R=X*J-J*X'-target;
    err=norm(R,'fro');
    if err<=tol
        solved=true;
        break
    end
    B=I+X;
    C=I-K;
    Jac=zeros(nd);
    for p=1:nd
        E=zeros(n);
        E(ii(p),jj(p))=1;
        E(jj(p),ii(p))=-1;
        DX=C\(E*B);
        DR=DX*J-J*DX';
        Jac(:,p)=DR(idx);
    end
    if rcond(Jac)<=1e-14
        error('um_moser_veselov_solve:NoSolution', ...
            'Implicit MV Jacobian is singular on selected Cayley branch.');
    end
    dx=-(Jac\R(idx));
    Dk=zeros(n);
    Dk(idx)=dx;
    Dk(sub2ind([n n],jj,ii))=-dx;
    improved=false;
    for ls=0:22
        Ktry=K+(2^(-ls))*Dk;
        Xtry=(I-Ktry)\(I+Ktry);
        Rtry=Xtry*J-J*Xtry'-target;
        if norm(Rtry,'fro')<(1-1e-5*2^(-ls))*err
            K=Ktry;
            improved=true;
            break
        end
    end
    if ~improved
        error('um_moser_veselov_solve:NoSolution', ...
            'No decrease on the locally selected Moser-Veselov branch.');
    end
end
if ~solved
    error('um_moser_veselov_solve:NoSolution', ...
        'No admissible converged SO(n) increment (try a smaller h).');
end
F=X';
if norm(F'*F-I,'fro')>2e-10*n ...
    || det(F)<=0 || norm(F'*J-J*F-target,'fro')>10*tol
    error('um_moser_veselov_solve:NoSolution', ...
        'Computed increment failed orthogonality or implicit residual.');
end
end
