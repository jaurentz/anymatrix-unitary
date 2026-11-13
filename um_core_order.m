function [order,steps]=um_core_order(n,representation,twists)
%UM_CORE_ORDER Adjacent Verblunsky-core order; 0 means terminal phase.
if nargin<3,twists=[];end
if nargin<2 || isempty(representation),representation='cmv';end
switch lower(char(representation))
    case 'hessenberg'
        order=1:n-1;steps=[order 0];
    case {'inverse_hessenberg','inverse'}
        order=n-1:-1:1;steps=[0 order];
    case 'cmv'
        odds=1:2:n-1;evens=2:2:n-1;
        order=[odds evens];
        % L*M: for odd n D is in L, preceding the even cores.
        % For even n D is in M and follows its even cores.
        if mod(n,2),steps=[odds 0 evens];
        else,steps=[odds evens 0];end
    case 'twisted'
        if isempty(twists),twists=mod((1:max(n-2,0))',2)==1;end
        validateattributes(twists,{'numeric','logical'},{'vector'});
        if numel(twists)~=max(0,n-2) || any(twists(:)~=0 & twists(:)~=1)
            error('js_core_product:Twists','TWISTS must have n-2 binary entries.');
        end
        order=1:min(n-1,1);
        for k=2:n-1
            if twists(k-1),order=[k order];else,order=[order k];end %#ok<AGROW>
        end
        steps=[order 0];
    otherwise
        error('js_core_product:Representation','Unknown representation %s.',char(representation));
end
end
