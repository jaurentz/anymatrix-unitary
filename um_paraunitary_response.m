function E=um_paraunitary_response(z,directions,U0)
%UM_PARAUNITARY_RESPONSE Evaluate the factored FIR response E(z).
d=size(directions,1);E=eye(d);
for j=1:size(directions,2)
    v=directions(:,j);
    E=E+(z^(-1)-1)*v*(v'*E);
end
E=U0*E;
end
