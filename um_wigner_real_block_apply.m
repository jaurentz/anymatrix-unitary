function Y=um_wigner_real_block_apply(X,flag,C,D)
%UM_WIGNER_REAL_BLOCK_APPLY Matrix-free real-harmonic rotation block.
%   Q=C'*Dcomplex*C is the real Condon--Shortley SO(3) representation,
%   where C is sparse unitary and Dcomplex is a Wigner operator handle.
%   Its mathematical transpose equals its adjoint because Q is real.
%   Complex right-hand sides remain complex-linear. Dropping small
%   imaginary roundoff is permitted ONLY for real input data.
if strcmpi(flag,'size'),Y=[size(C,1) size(C,1)];return;end
flag=lower(char(flag));
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_wigner_real_block_apply:Flag', ...
        'Use notransp, transp, ctransp or size.');
end
if strcmp(flag,'notransp')
    mode='notransp';
else
    mode='ctransp';
end
Y=C'*D(C*X,mode);
if isreal(X),Y=real(Y);end
end
