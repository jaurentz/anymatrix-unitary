function Y=um_dct_apply(X,flag,n,kind)
%UM_DCT_APPLY FFT-based DCT-I through IV; forward and real transpose.
%   Implements orthonormal transforms without forming cosine matrices.
%   For complex RHS, cosine transforms remain complex linear.
if strcmpi(flag,'size'),Y=[n n];return;end
switch lower(char(flag))
    case 'notransp'
    case {'transp','ctransp'}
        if kind==2,kind=3;elseif kind==3,kind=2;end
    otherwise
        error('um_dct_apply:Flag','Unknown transform action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_dct_apply:Dimension','Input must have n rows.');
end
if n==1
    % DCT-II, III and IV all reduce to the scalar identity.
    Y=X;
    return
end
wasReal=isreal(X);
switch kind
    case 1
        c=ones(n,1);c([1 n])=1/sqrt(2);
        V=c.*X;
        V([1 n],:)=2*V([1 n],:);
        W=fft([V;flipud(V(2:end-1,:))],[],1);
        Y=(c/sqrt(2*(n-1))).*W(1:n,:);
    case 2
        c=ones(n,1);c(1)=1/sqrt(2);
        k=(0:n-1)';
        W=fft([X;flipud(X)],[],1);
        Y=(c.*exp(-1i*pi*k/(2*n))/sqrt(2*n)).*W(1:n,:);
    case 3
        % Full complex-linear inverse of the DCT-II FFT relation.
        c=ones(n,1);c(1)=1/sqrt(2);
        k=(0:n-1)';
        W=(sqrt(2*n)*exp(1i*pi*k/(2*n))./c).*X;
        tail=exp(-1i*pi*flipud(k(2:end))/n).*flipud(W(2:end,:));
        fullSpectrum=[W;zeros(1,size(X,2),'like',X);tail];
        T=ifft(fullSpectrum,[],1);
        Y=T(1:n,:);
    case 4
        W=um_dct_apply([X;zeros(n,size(X,2),'like',X)], ...
            'notransp',2*n,2);
        Y=sqrt(2)*W(2:2:2*n,:);
end
if wasReal,Y=real(Y);end
end
