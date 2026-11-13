function Y=um_dst_apply(X,flag,n,kind)
%UM_DST_APPLY Orthonormal DST-I through IV via FFT and shared DCT kernels.
if strcmpi(flag,'size'),Y=[n n];return;end
switch lower(char(flag))
    case 'notransp'
    case {'transp','ctransp'}
        if kind==2,kind=3;elseif kind==3,kind=2;end
    otherwise
        error('um_dst_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_dst_apply:Dimension','Input must have n rows.');
end
s=(-1).^(0:n-1)';
switch kind
    case 1
        zero=zeros(1,size(X,2),'like',X);
        W=fft([zero;X;zero;-flipud(X)],[],1);
        Y=(1i/sqrt(2*(n+1)))*W(2:n+1,:);
        if isreal(X),Y=real(Y);end
    case 2
        Y=flipud(um_dct_apply(s.*X,'notransp',n,2));
    case 3
        Y=s.*um_dct_apply(flipud(X),'transp',n,2);
    case 4
        Y=flipud(um_dct_apply(s.*X,'notransp',n,4));
end
end
