function v=um_finite_koopman_cycle_moments(lengths,N,k)
%UM_FINITE_KOOPMAN_CYCLE_MOMENTS Tr(U^k)/N as fixed-point fraction.
validateattributes(k,{'numeric'},{'real','finite','integer'});
v=zeros(size(k));
for j=1:numel(k)
    v(j)=sum(lengths(mod(k(j),lengths)==0))/N;
end
end
