function tf=um_is_operator(mode)
%UM_IS_OPERATOR Return true for an implicit operator output request.
if nargin==0 || isempty(mode),tf=false;return;end
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('um_is_operator:Mode','Output mode must be matrix or operator.');
end
switch lower(char(mode))
    case {'matrix','explicit'},tf=false;
    case {'operator','implicit'},tf=true;
    otherwise,error('um_is_operator:Mode','Unknown output mode %s.',char(mode));
end
end
