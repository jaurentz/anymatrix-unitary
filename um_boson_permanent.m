function v=um_boson_permanent(A)
%UM_BOSON_PERMANENT Small permanent by Ryser's inclusion-exclusion.
p=size(A,1);
if p==0,v=1;return;end
v=0;
for mask=1:(2^p-1)
    cols=logical(bitget(mask,1:p));
    s=sum(A(:,cols),2);
    v=v+(-1)^(p-sum(cols))*prod(s);
end
end
