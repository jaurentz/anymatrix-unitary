function Y=um_wavelet_apply(X,flag,n,levels)
%UM_WAVELET_APPLY Haar analysis transform and its conjugate transpose.
if strcmpi(flag,'size'),Y=[n n];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_wavelet_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_wavelet_apply:Dimension','Input must have n rows.');
end
Y=X;
if adj,levelsToApply=levels:-1:1;else,levelsToApply=1:levels;end
for ell=levelsToApply
    m=n/2^(ell-1);h=m/2;
    T=Y(1:m,:);
    if adj
        low=T(1:h,:);high=T(h+1:end,:);
        Y(1:2:m,:)=(low+high)/sqrt(2);
        Y(2:2:m,:)=(low-high)/sqrt(2);
    else
        a=T(1:2:m,:);b=T(2:2:m,:);
        Y(1:m,:)=[(a+b);(a-b)]/sqrt(2);
    end
end
end
