function m=um_optical_mesh_moments(k,z)
%UM_OPTICAL_MESH_MOMENTS Exact trace moments for diagonal mesh.
validateattributes(k,{'numeric'},{'real','integer','finite'});
shape=size(k);
if isempty(k),m=zeros(shape);return;end
m=reshape(mean(z(:).^reshape(k,1,[]),1),shape);
end
