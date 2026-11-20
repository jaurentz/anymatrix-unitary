function m=um_skew_moments(z,k)
%UM_SKEW_MOMENTS Normalized exact algebraic constant-flow trace moments.
if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
        ||any(k(:)~=round(k(:)))
    error('um_skew_moments:Moments','K must contain finite integers.');
end
m=zeros(size(k));
for j=1:numel(k),m(j)=mean(z.^k(j));end
end
