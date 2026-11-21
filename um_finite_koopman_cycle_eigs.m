function z=um_finite_koopman_cycle_eigs(lengths)
%UM_FINITE_KOOPMAN_CYCLE_EIGS Roots of unity from each orbit.
N=sum(lengths);z=complex(zeros(N,1));start=1;
for j=1:numel(lengths)
    ell=lengths(j);
    z(start:start+ell-1)=exp(2i*pi*(0:ell-1)'/ell);
    start=start+ell;
end
end
