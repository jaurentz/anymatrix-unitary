function Y=um_skew_midpoint_apply(X,flag,n,stages,method)
%UM_SKEW_MIDPOINT_APPLY Factor-wise orthogonal forward and adjoint.
if ~(ischar(flag)||(isstring(flag)&&isscalar(flag)))
    error('um_skew_midpoint_apply:Flag','Invalid operator flag.');
end
flag=lower(char(flag));
if strcmp(flag,'size'),Y=[n n];return;end
if ~isnumeric(X)||ndims(X)>2||size(X,1)~=n
    error('um_skew_midpoint_apply:Dimension','Expected N input rows.');
end
if ~any(strcmp(flag,{'notransp','transp','ctransp'}))
    error('um_skew_midpoint_apply:Flag','Unknown action.');
end
if issparse(X),X=full(X);end
Y=X;
if strcmp(flag,'notransp')
    range=1:numel(stages);
else
    range=numel(stages):-1:1;
end
for j=range
    if strcmp(method,'exp_midpoint')
        Q=stages{j};
        if strcmp(flag,'notransp'),Y=Q*Y;
        else,Y=Q'*Y;end
    else
        % (I-C)^(-1)(I+C)=2(I-C)^(-1)-I, since C is skew.
        % Its adjoint is 2(I+C)^(-1)-I, reverse factor order.
        previous=Y;
        if strcmp(flag,'notransp')
            Y=2*(stages{j}.minus\Y)-previous;
        else
            Y=2*(stages{j}.plus\Y)-previous;
        end
    end
end
end
