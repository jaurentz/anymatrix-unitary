function out=um_exterior_spectral(z,p,mode,k)
%UM_EXTERIOR_SPECTRAL Analytic subset-product phases and trace moments.
d=numel(z);N=nchoosek(d,p);
switch lower(char(mode))
    case 'eigenvalues'
        if N>2e6
            error('um_exterior_spectral:Materialization', ...
                'Too many eigenvalues; use exact finite.moments.');
        end
        if p==0,out=1;return;end
        if p==d,out=prod(z);return;end
        idx=nchoosek(1:d,p);
        out=prod(z(idx),2);
        out=out(:);
    case 'moments'
        if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
                ||any(k(:)~=round(k(:)))
            error('um_exterior_spectral:Moments','K must be finite integers.');
        end
        out=zeros(size(k));
        for t=1:numel(k)
            c=zeros(p+1,1);c(1)=1;
            for j=1:d
                w=z(j)^k(t);
                for r=min(p,j):-1:1
                    c(r+1)=c(r+1)+w*c(r);
                end
            end
            out(t)=c(p+1)/N;
        end
    otherwise
        error('um_exterior_spectral:Mode','Unknown spectral mode.');
end
end
