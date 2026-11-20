function Y=um_graph_connection_cayley(X,flag,N,steps,plusFactor,minusFactor)
%UM_GRAPH_CONNECTION_CAYLEY Unitary rational approximation to exp(-itL).
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_graph_connection_cayley:Flag','Unknown action flag.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[N N];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=N
    error('um_graph_connection_cayley:Dimension','Expected N*Q rows.');
end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_graph_connection_cayley:Flag','Unknown action.');
end
if issparse(X),X=full(X);end
Y=X;
for j=1:steps
    old=Y;
    if strcmp(flag,'notransp')
        Y=2*(plusFactor\Y)-old;
    else
        Y=2*(minusFactor\Y)-old;
    end
end
end
