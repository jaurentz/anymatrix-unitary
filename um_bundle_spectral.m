function out=um_bundle_spectral(permutation,links,mode,k)
%UM_BUNDLE_SPECTRAL Exact finite phases and trace moments by cycle.
n=numel(permutation);d=size(links,1);N=n*d;
seen=false(n,1);lens=[];H={};
for start=1:n
    if seen(start),continue;end
    ell=0;j=start;product=eye(d);
    while ~seen(j)
        seen(j)=true;ell=ell+1;
        product=links(:,:,j)*product;
        j=permutation(j);
    end
    lens(end+1)=ell; %#ok<AGROW>
    H{end+1}=product; %#ok<AGROW>
end
switch lower(char(mode))
    case 'eigenvalues'
        if N>2e6
            error('um_bundle_spectral:Materialization','Use moments for large N.');
        end
        out=complex(zeros(N,1));used=0;
        for c=1:numel(lens)
            ell=lens(c);z=eig(H{c});
            for j=1:d
                out(used+(1:ell))=exp(1i*(angle(z(j))+2*pi*(0:ell-1)')/ell);
                used=used+ell;
            end
        end
    case 'moments'
        if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
                ||any(k(:)~=round(k(:)))
            error('um_bundle_spectral:Moments','K must be finite integers.');
        end
        out=zeros(size(k));
        for r=1:numel(k)
            t=0;
            for c=1:numel(lens)
                ell=lens(c);
                if mod(k(r),ell)~=0,continue;end
                h=k(r)/ell;
                if h>=0,W=H{c}^h;else,W=(H{c}')^(-h);end
                t=t+ell*trace(W);
            end
            out(r)=t/N;
        end
    otherwise
        error('um_bundle_spectral:Mode','Invalid spectral mode.');
end
end
