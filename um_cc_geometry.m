function [inputs,outputs,gates]=um_cc_geometry(nx,ny,theta)
%UM_CC_GEOMETRY Disjoint two-channel ports of the directed CC network.
%   INPUTS and OUTPUTS are 2-by-(2*nx*ny), with disjoint column pairs.
%   GATES is 4-by-(2*nx*ny): each column is a column-major 2x2 gate.
%   The same geometry serves sparse construction and implicit actions.
a=[cos(theta(1)),sin(theta(1));-sin(theta(1)),cos(theta(1))];
b=[sin(theta(2)),cos(theta(2));-cos(theta(2)),sin(theta(2))];
X=2*nx;Y=2*ny;
id=@(x,y)1+mod(x,X)+X*mod(y,Y);
m=2*nx*ny;
inputs=zeros(2,m);
outputs=zeros(2,m);
gates=zeros(4,m);
q=0;
for j=0:nx-1
    for k=0:ny-1
        q=q+1;
        inputs(:,q)=[id(2*j,2*k);id(2*j+1,2*k-1)];
        outputs(:,q)=[id(2*j+1,2*k);id(2*j,2*k-1)];
        gates(:,q)=a(:);
        q=q+1;
        inputs(:,q)=[id(2*j+2,2*k+1);id(2*j+1,2*k)];
        outputs(:,q)=[id(2*j+1,2*k+1);id(2*j+2,2*k)];
        gates(:,q)=b(:);
    end
end
end
