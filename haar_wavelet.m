function [U,properties] = haar_wavelet(n,levels,output)
%HAAR_WAVELET Finite orthonormal Haar wavelet analysis transform.
%   U=HAAR_WAVELET(N,LEVELS) forms a real orthogonal N-by-N transform.
%   N is a power of two. Each of LEVELS decomposition steps computes
%   averages (x_(2j-1)+x_(2j))/sqrt(2) and details with minus sign,
%   operating on the averages from the previous step. Default LEVELS
%   is log2(N); LEVELS=0 returns I. Coarse coefficients appear first.
%   Reference: I. Daubechies, Ten Lectures on Wavelets, SIAM, 1992.
%   https://doi.org/10.1137/1.9781611970104
properties={'square','unitary','normal','orthogonal','real','sparse', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
maxLevels=round(log2(n));
if 2^maxLevels~=n,error('haar_wavelet:Order','N must be a power of two.');end
if nargin<2 || isempty(levels),levels=maxLevels;end
validateattributes(levels,{'numeric'},{'scalar','integer','>=',0,'<=',maxLevels});
if nargin<3,output='matrix';end
if um_is_operator(output)
    U=@(X,flag)um_wavelet_apply(X,flag,n,levels);
    return
end
U=speye(n);
for ell=1:levels
    m=n/2^(ell-1);h=m/2;
    ii=[1:h,1:h,h+(1:h),h+(1:h)];
    jj=[1:2:m,2:2:m,1:2:m,2:2:m];
    vv=[ones(1,3*h),-ones(1,h)]/sqrt(2);
    S=sparse(ii,jj,vv,m,m);
    U(1:m,:)=S*U(1:m,:);
end
end
