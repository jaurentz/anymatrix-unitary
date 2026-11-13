function Y=um_hartley_apply(X,flag,n)
%UM_HARTLEY_APPLY FFT action of the symmetric Hartley involution.
if strcmpi(flag,'size'),Y=[n n];return;end
if ~any(strcmpi(flag,{'notransp','transp','ctransp'}))
    error('um_hartley_apply:Flag','Unknown transform action.');
end
if ~isnumeric(X) || ndims(X)>2 || size(X,1)~=n
    error('um_hartley_apply:Dimension','Input must have n rows.');
end
% H=[(1+i)F+(1-i)F*]/2, with F the normalized DFT.
Y=((1+1i)*um_fourier_apply(X,'notransp',0)+ ...
   (1-1i)*um_fourier_apply(X,'transp',0))/2;
if isreal(X),Y=real(Y);end
end
