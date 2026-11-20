function out=um_paraunitary_spectral(L,directions,U0,mode,k)
%UM_PARAUNITARY_SPECTRAL Exact periodic Fourier-block phase reduction.
d=size(directions,1);n=d*L;
if strcmpi(mode,'eigenvalues')&&n>2e6
    error('um_paraunitary_spectral:Materialization', ...
        'Use Fourier moments for orders above 2e6.');
end
if strcmpi(mode,'moments')
    if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
            ||any(k(:)~=round(k(:)))
        error('um_paraunitary_spectral:Moments', ...
            'Fourier moment orders must be finite integers.');
    end
end
if strcmpi(mode,'eigenvalues')
    out=complex(zeros(n,1));
elseif strcmpi(mode,'moments')
    out=zeros(size(k));out=complex(out);
else
    error('um_paraunitary_spectral:Mode','Unknown spectral mode.');
end
for j=0:L-1
    E=um_paraunitary_response(exp(2i*pi*j/L),directions,U0);
    z=eig(E);
    z=z./abs(z);
    if strcmpi(mode,'eigenvalues')
        out(j*d+(1:d))=z;
    else
        for t=1:numel(k)
            out(t)=out(t)+sum(z.^k(t))/n;
        end
    end
end
end
