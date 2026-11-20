function M=um_block_cmv_moments(nodes,weights,k,mode)
%UM_BLOCK_CMV_MOMENTS Matrix-valued and scalar trace Fourier moments.
if ~isnumeric(k)||~isreal(k)||any(~isfinite(k(:))) ...
        ||any(k(:)~=round(k(:)))
    error('um_block_cmv_moments:Order','Fourier orders must be finite integers.');
end
if strcmp(mode,'trace')
    M=zeros(size(k));
    for j=1:numel(k)
        M(j)=mean(nodes.^k(j));
    end
elseif strcmp(mode,'cyclic')
    q=size(weights,1);
    M=complex(zeros(q,q,numel(k)));
    for j=1:numel(k)
        Z=zeros(q);
        for t=1:numel(nodes)
            Z=Z+nodes(t)^k(j)*weights(:,:,t);
        end
        M(:,:,j)=Z;
    end
else
    error('um_block_cmv_moments:Mode','Invalid moments mode.');
end
end
