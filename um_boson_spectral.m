function out=um_boson_spectral(z,p,mode,k)
%UM_BOSON_SPECTRAL Occupation products and complete homogeneous moments.
d=numel(z);N=nchoosek(d+p-1,p);
switch lower(char(mode))
    case 'eigenvalues'
        if N>2e6
            error('um_boson_spectral:Materialization', ...
                'Too many eigenvalues; use finite.moments instead.');
        end
        if p==0,out=1;return;end
        if d==1,out=z^p;return;end
        occ=um_boson_occupations(d,p);
        out=prod(bsxfun(@power,z(:).',occ),2);
        out=out(:);
    case 'moments'
        if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
                ||any(k(:)~=round(k(:)))
            error('um_boson_spectral:Moments','K must contain finite integers.');
        end
        out=zeros(size(k));
        for t=1:numel(k)
            h=zeros(p+1,1);h(1)=1;
            for j=1:d
                w=z(j)^k(t);
                for r=1:p
                    % Ascending degree admits repeated occupation of mode j.
                    h(r+1)=h(r+1)+w*h(r);
                end
            end
            out(t)=h(p+1)/N;
        end
    otherwise
        error('um_boson_spectral:Mode','Invalid spectral mode.');
end
end
