function Y=um_spin_sampling_apply(X,flag,data)
%UM_SPIN_SAMPLING_APPLY Exact bandlimited spin-s sphere synthesis/adjoint.
%   DATA has fields L, M=2L-1, Ncoef=L^2-s^2, B_m
%   (L-by-[L-max(|s|,|m|)] weighted latitude blocks) and
%   global columns selected by fixed m. Both actions are FFT
%   in longitude and dense Gauss latitude recurrence blocks.
%   The rectangular sampling isometry maps Ncoef coefficients
%   to L*(2L-1) weighted samples. The size interface exposes
%   this RECTANGULAR shape. Storage/setup O(L^3), apply
%   O(L^3+L^2 log L) per RHS, not a square unitary.
L=data.L;M=data.M;C=data.Ncoef;rows=L*M;
if strcmpi(flag,'size'),Y=[rows C];return;end
if ~isnumeric(X)||ndims(X)>2
    error('um_spin_sampling_apply:Dimension','Numeric matrix input required.');
end
switch lower(char(flag))
    case 'notransp'
        if size(X,1)~=C
            error('um_spin_sampling_apply:Dimension', ...
                'Synthesis requires L^2-s^2 coefficient rows.');
        end
        Y=complex(zeros(rows,size(X,2)));
        for rhs=1:size(X,2)
            freq=complex(zeros(L,M));
            for m=-(L-1):(L-1)
                j=m+L;ids=data.cols{j};
                if isempty(ids),continue;end
                freq(:,mod(m,M)+1)=data.blocks{j}*X(ids,rhs);
            end
            samples=sqrt(M)*ifft(freq,[],2);
            Y(:,rhs)=samples(:);
        end
    case {'transp','ctransp'}
        if size(X,1)~=rows
            error('um_spin_sampling_apply:Dimension', ...
                'Analysis requires L*(2L-1) sample rows.');
        end
        Y=complex(zeros(C,size(X,2)));
        for rhs=1:size(X,2)
            samples=reshape(X(:,rhs),L,M);
            freq=fft(samples,[],2)/sqrt(M);
            for m=-(L-1):(L-1)
                j=m+L;ids=data.cols{j};
                if isempty(ids),continue;end
                Y(ids,rhs)=data.blocks{j}'*freq(:,mod(m,M)+1);
            end
        end
    otherwise
        error('um_spin_sampling_apply:Flag', ...
            'Use notransp, transp, ctransp, or size.');
end
end
