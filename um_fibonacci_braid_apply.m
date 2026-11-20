function Y=um_fibonacci_braid_apply(X,flag,d,sig,word)
%UM_FIBONACCI_BRAID_APPLY Sparse factors for signed braid words.
%   Each positive word index applies sigma_i. Negative applies sigma_i'.
%   'transp' follows the collection convention for the real/complex
%   adjoint ('ctransp'); ordinary nonconjugating transpose is not used.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_fibonacci_braid_apply:Flag','Invalid action flag.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[d d];return;end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_fibonacci_braid_apply:Flag','Unknown action.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=d
    error('um_fibonacci_braid_apply:Dimension','Expected fusion-space rows.');
end
Y=X;
if strcmp(flag,'notransp')
    indices=1:numel(word);
    adjoint=false;
else
    indices=numel(word):-1:1;
    adjoint=true;
end
for k=indices
    j=word(k);
    G=sig{abs(j)};
    if xor(j<0,adjoint),G=G';end
    Y=G*Y;
end
end
