function C=um_graph_coin_matrix(A,src,dst,starts,coins,weighted)
%UM_GRAPH_COIN_MATRIX Block-local unitary coins on ordered outgoing arcs.
%   COINS: 'grover', 'fourier', 'identity', or cell of local unitaries.
%   WEIGHTED true applies edge weights to the Grover reflection vector.
n=numel(starts)-1;m=numel(src);
if iscell(coins)
    if numel(coins)~=n
        error('um_graph_coin_matrix:Coins','Need one coin per vertex.');
    end
else
    if ~(ischar(coins)||isstring(coins)) || ~isscalar(string(coins))
        error('um_graph_coin_matrix:Coins','Invalid coin type.');
    end
    coins=lower(char(coins));
    if ~any(strcmp(coins,{'grover','fourier','identity'}))
        error('um_graph_coin_matrix:Coins','Unknown coin type.');
    end
end
C=spalloc(m,m,sum(diff(starts).^2));
for v=1:n
    ix=starts(v):starts(v+1)-1;d=numel(ix);
    if iscell(coins)
        Q=coins{v};
        if ~(isnumeric(Q)&&isequal(size(Q),[d d])&&all(isfinite(Q(:))) ...
                &&norm(Q'*Q-eye(d),'fro')<=1e-10*max(1,d))
            error('um_graph_coin_matrix:Coins', ...
                  'Every supplied local coin must be a unitary d-by-d matrix.');
        end
    else
        switch coins
            case 'identity'
                Q=eye(d);
            case 'fourier'
                Q=fft(eye(d))/sqrt(d);
            otherwise
                if weighted
                    w=full(A(sub2ind(size(A),src(ix),dst(ix))));
                else
                    w=ones(d,1);
                end
                s=sqrt(w(:)/sum(w));
                Q=2*(s*s')-eye(d);
        end
    end
    C(ix,ix)=Q;
end
end
