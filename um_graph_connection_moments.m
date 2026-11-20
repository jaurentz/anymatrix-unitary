function out=um_graph_connection_moments(z,k)
%UM_GRAPH_CONNECTION_MOMENTS Uniform Fourier moments of twisted cycle.
if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
        ||any(k(:)~=round(k(:)))
    error('um_graph_connection_moments:Order', ...
        'Moment indices must be finite integers.');
end
out=zeros(size(k));
for j=1:numel(k),out(j)=mean(z.^k(j));end
end
