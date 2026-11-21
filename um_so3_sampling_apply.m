function Y=um_so3_sampling_apply(X,flag,data)
%UM_SO3_SAMPLING_APPLY FFT-separated Peter--Weyl synthesis/analysis.
%   Gauss-Legendre beta sampling and (2L-1)^2 uniform alpha/gamma.
%   If M=2L-1, rows L*M^2, columns sum_(ell<L)(2ell+1)^2.
%   Rows beta FASTEST, then alpha, then gamma; columns ell,m,n.
%   A'*A=I in exact arithmetic, while A is rectangular except L=1.
%   For fixed (m,n), weighted Wigner latitude block B_(m,n) has
%   dimension L-by-(L-max(abs(m),abs(n))) and is built by
%   UM_WIGNER_JACOBI. The remaining two angles are orthonormal FFTs.
%   Storage O(L^4), per RHS cost O(L^4+L^3 log L), compared with
%   O(L^6) for explicitly storing/applying dense sampling A.
%
%   Haar measure normalized as d(alpha) sin(beta)d(beta)d(gamma)/(8pi^2),
%   basis functions sqrt(2ell+1) D^ell_(m,n)(alpha,beta,gamma),
%   D=exp(-im alpha)*d^ell_(m,n)(beta)*exp(-in gamma).
%   This is NOT the distinct Wigner rotation representation
%   D^ell(R) acting on a (2ell+1)-coefficient vector.
L=data.L;M=data.M;C=data.Ncoef;rows=L*M*M;
if strcmpi(flag,'size'),Y=[rows C];return;end
if ~isnumeric(X)||ndims(X)>2
    error('um_so3_sampling_apply:Dimension','Numeric matrix input required.');
end
switch lower(char(flag))
    case 'notransp'
        if size(X,1)~=C
            error('um_so3_sampling_apply:Dimension', ...
                'Synthesis requires the full bandlimited SO(3) coefficient count.');
        end
        Y=complex(zeros(rows,size(X,2)));
        for rhs=1:size(X,2)
            freq=complex(zeros(L,M,M));
            for m=-(L-1):L-1
                i=mod(m,M)+1;
                for n=-(L-1):L-1
                    j=mod(n,M)+1;
                    ids=data.cols{i,j};
                    if isempty(ids),continue;end
                    freq(:,i,j)=data.blocks{i,j}*X(ids,rhs);
                end
            end
            samples=fft(fft(freq,[],2),[],3)/M;
            Y(:,rhs)=samples(:);
        end
    case {'transp','ctransp'}
        if size(X,1)~=rows
            error('um_so3_sampling_apply:Dimension', ...
                'Analysis requires L*(2L-1)^2 sample rows.');
        end
        Y=complex(zeros(C,size(X,2)));
        for rhs=1:size(X,2)
            samples=reshape(X(:,rhs),L,M,M);
            freq=M*ifft(ifft(samples,[],2),[],3);
            for m=-(L-1):L-1
                i=mod(m,M)+1;
                for n=-(L-1):L-1
                    j=mod(n,M)+1;
                    ids=data.cols{i,j};
                    if isempty(ids),continue;end
                    Y(ids,rhs)=data.blocks{i,j}'*freq(:,i,j);
                end
            end
        end
    otherwise
        error('um_so3_sampling_apply:Flag', ...
            'Use notransp, transp, ctransp, or size.');
end
end
