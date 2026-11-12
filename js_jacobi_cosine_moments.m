function moments = js_jacobi_cosine_moments(a,b,kmax)
%JS_JACOBI_COSINE_MOMENTS Exact finite cosine moments of Jacobi Szego lift.
%   MOMENTS(k+1)=integral cos(k*theta) dmu(theta).
%   Uses E[T_k(t)] with t=2*Y-1, Y~Beta(b+1,a+1). Polynomial arithmetic
%   is intended for modest k (e.g., up to 8), not high-order stable DOS.
validateattributes(kmax,{'numeric'},{'scalar','integer','nonnegative'});
betaMom=ones(kmax+1,1);
for k=1:kmax
    betaMom(k+1)=betaMom(k)*(b+k)/(a+b+k+1);
end
moments=zeros(kmax+1,1);
t0=1; t1=[-1,2];
for k=0:kmax
    if k==0
        pk=t0;
    elseif k==1
        pk=t1;
    else
        temp=[0,4*t1]-[2*t1,0]; % 2*(2Y-1)*T_(k-1)
        pk=temp-[t0,zeros(1,numel(temp)-numel(t0))];
        t0=t1;t1=pk;
    end
    moments(k+1)=pk*betaMom(1:numel(pk));
end
end
