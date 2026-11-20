function C=um_paraunitary_taps(directions,U0)
%UM_PARAUNITARY_TAPS Laurent coefficient matrices of E(z).
d=size(directions,1);m=size(directions,2);
C=zeros(d,d,m+1);
if ~isreal(directions)||~isreal(U0),C=complex(C);end
C(:,:,1)=eye(d);
for j=1:m
    v=directions(:,j);P=v*v';
    Q=eye(d)-P;
    next=zeros(d,d,m+1,'like',C);
    for r=0:j
        if r<j
            next(:,:,r+1)=next(:,:,r+1)+Q*C(:,:,r+1);
        end
        if r>0
            next(:,:,r+1)=next(:,:,r+1)+P*C(:,:,r);
        end
    end
    C=next;
end
for r=1:m+1,C(:,:,r)=U0*C(:,:,r);end
end
