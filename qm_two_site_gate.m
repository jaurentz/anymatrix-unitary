function A=qm_two_site_gate(n,gate,site1,site2)
%QM_TWO_SITE_GATE Embed a two-qubit gate at ordered sites in n qubits.
%   Computational basis: site 1 is the most-significant tensor factor.
%   The 4-by-4 local basis is |00>, |01>, |10>, |11>.
%   Works also for a periodic wrap bond with ordered sites (n,1).
%   Internal helper for Floquet generators (not an Anymatrix generator).
validateattributes(n,{'numeric'},{'scalar','integer','>=',2,'<=',20,'finite'});
validateattributes(gate,{'numeric'},{'size',[4,4],'finite'});
validateattributes(site1,{'numeric'},{'scalar','integer','>=',1,'<=',n});
validateattributes(site2,{'numeric'},{'scalar','integer','>=',1,'<=',n});
if site1==site2
    error('qm_two_site_gate:Sites','Sites must be distinct.');
end
N=2^n;
row=zeros(4*N,1);
col=zeros(4*N,1);
val=complex(zeros(4*N,1));
t=0;
for state=0:N-1
    s=uint64(state);
    b1=double(bitget(s,n-site1+1));
    b2=double(bitget(s,n-site2+1));
    input=1+2*b1+b2;
    for output=0:3
        c1=floor(output/2);
        c2=mod(output,2);
        target=bitset(s,n-site1+1,c1);
        target=bitset(target,n-site2+1,c2);
        t=t+1;
        row(t)=double(target)+1;
        col(t)=state+1;
        val(t)=gate(output+1,input);
    end
end
A=sparse(row,col,val,N,N);
end
