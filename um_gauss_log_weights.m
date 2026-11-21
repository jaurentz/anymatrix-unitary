function logw=um_gauss_log_weights(nodes,diagonal,offdiagonal)
%UM_GAUSS_LOG_WEIGHTS Logarithmic Gaussian Christoffel weights.
%   LOGW(j)=-log(sum(k=0:n-1, p_k(nodes(j))^2)), where p_k
%   are the orthonormal polynomials of the positive-offdiagonal
%   Jacobi recurrence p_0=1, p_-1=0. Thus exp(LOGW) equals
%   NORMALIZED Gauss quadrature weights whenever they are
%   representable. Unlike squared first eigenvector components,
%   these log-weights can resolve values vastly below realmin.
%
%   Each forward recurrence is rescaled dynamically before overflow
%   and squared sums are accumulated using log-sum-exp. O(n^2)
%   arithmetic, O(n) output plus scalar work. The log sums are
%   approximate in floating point; no rigorous relative accuracy
%   is claimed for arbitrarily large order, clustered quadrature
%   nodes or near-singular Jacobi parameters.
x=nodes(:);d=diagonal(:);e=offdiagonal(:);
n=numel(d);logw=zeros(numel(x),1);
for j=1:numel(x)
    previous=0;
    current=1;
    logScale=0;
    logNorm2=0; % includes p_0=1, hence exact initial sum=1
    for k=1:n-1
        if k==1,minus=0;else,minus=e(k-1)*previous;end
        next=((x(j)-d(k))*current-minus)/e(k);
        if ~isfinite(next)
            error('um_gauss_log_weights:Recurrence', ...
                'Non-finite Jacobi recurrence at node %d.',j);
        end
        if max(abs(next),abs(current))>1e90
            previous=previous*1e-90;
            current=current*1e-90;
            next=next*1e-90;
            logScale=logScale+90*log(10);
        end
        if next~=0
            logterm=2*(log(abs(next))+logScale);
            m=max(logNorm2,logterm);
            logNorm2=m+log(exp(logNorm2-m)+exp(logterm-m));
        end
        previous=current;
        current=next;
    end
    logw(j)=-logNorm2;
end
end
