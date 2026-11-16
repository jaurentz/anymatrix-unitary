function Y=um_graph_walk_apply(X,flag,C,reverse,phase,kind)
%UM_GRAPH_WALK_APPLY Matrix-free coined and metric-graph arc actions.
%   COINED U=S*C; METRIC U=diag(phase)*C*S, S the arc reversal.
m=numel(reverse);
if strcmpi(flag,'size'),Y=[m m];return;end
adj=any(strcmpi(flag,{'transp','ctransp'}));
if ~adj && ~strcmpi(flag,'notransp')
    error('um_graph_walk_apply:Flag','Unknown operator action.');
end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=m
    error('um_graph_walk_apply:Dimension','Input has wrong row count.');
end
if strcmp(kind,'coined')
    if adj,Y=C'*X(reverse,:);
    else,Y=C*X;Y=Y(reverse,:);end
elseif strcmp(kind,'metric')
    if adj
        Y=C'*(conj(phase(:)).*X);
        Y=Y(reverse,:);
    else
        Y=phase(:).*(C*X(reverse,:));
    end
else
    error('um_graph_walk_apply:Kind','Unknown walk kind.');
end
end
