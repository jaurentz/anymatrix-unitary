function [cycles,lengths]=um_finite_koopman_cycles(perm)
%UM_FINITE_KOOPMAN_CYCLES Canonical disjoint-cycle decomposition.
N=numel(perm);
visited=false(N,1);cycles=cell(0,1);lengths=zeros(0,1);
for j=1:N
    if visited(j),continue;end
    k=j;orbit=zeros(0,1);
    while ~visited(k)
        visited(k)=true;
        orbit(end+1,1)=k; %#ok<AGROW>
        k=perm(k);
    end
    if k~=j
        error('um_finite_koopman_cycles:Permutation', ...
            'Input is not a bijection.');
    end
    cycles{end+1,1}=orbit; %#ok<AGROW>
    lengths(end+1,1)=numel(orbit); %#ok<AGROW>
end
end
