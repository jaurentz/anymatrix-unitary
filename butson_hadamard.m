function [U,properties]=butson_hadamard(n,construction,output)
%BUTSON_HADAMARD Algebraic Butson matrices with exact phase tables.
%   U=BUTSON_HADAMARD(N,CONSTRUCTION) returns a normalized unitary
%   complex Hadamard H/sqrt(N), where H's entries are q-th roots of
%   unity (H is BH(N,q)). 'tao' is the isolated 6-by-6 BH(6,3)
%   matrix S_6 of Tao, reported by Tadej and Zyczkowski.
%   'paley' is the Paley type-I real BH(p+1,2) for a prime
%   p=N-1 congruent to 3 mod 4 (e.g. orders 4,8,12,20,24).
%   For Paley, the p-by-p core Q(i,j) is the Legendre symbol
%   chi(i-j) for i~=j, zero on the diagonal; then
%       H=[1,ones(1,p);ones(p,1),Q-eye(p)].
%   In each case H*H'=N*I in exact arithmetic.
%   OUTPUT='exponents' returns integer E and root_order q with
%   H(j,k)=exp(2*pi*1i*E(j,k)/q), including exact signs for q=2.
%   'spectral' labels generic eigenphases unknown: equivalence
%   under row/column phases and permutations is NOT similarity.
%   No fast implicit action is claimed for these dense constructions.
%   [92] https://doi.org/10.1007/s11080-006-8220-2
%   [121] https://doi.org/10.1002/sapm1933121311
properties={'square','unitary','normal','complex','scalable','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(construction)
    if n==6,construction='tao';else,construction='paley';end
end
if ~(ischar(construction)||(isstring(construction)&&isscalar(construction)))
    error('butson_hadamard:Construction','Invalid construction.');
end
construction=lower(char(construction));
if n>512
    error('butson_hadamard:Order','Explicit Butson orders capped at N=512.');
end
switch construction
    case 'tao'
        if n~=6
            error('butson_hadamard:Order','Tao BH(6,3) requires N=6.');
        end
        q=3;
        E=[0 0 0 0 0 0; ...
           0 0 1 1 2 2; ...
           0 1 0 2 2 1; ...
           0 1 2 0 1 2; ...
           0 2 2 1 0 1; ...
           0 2 1 2 1 0];
    case 'paley'
        p=n-1;
        if p<3||mod(p,4)~=3||~isprime(p)
            error('butson_hadamard:Order', ...
                'Paley type I requires N=p+1 and p prime, p=3 mod 4.');
        end
        q=2;
        residues=false(p,1);
        for t=1:p-1
            residues(mod(t*t,p)+1)=true;
        end
        Q=zeros(p);
        for i=0:p-1
            for j=0:p-1
                if i==j,continue;end
                r=mod(i-j,p);
                Q(i+1,j+1)=2*double(residues(r+1))-1;
            end
        end
        H=[1,ones(1,p);ones(p,1),Q-eye(p)];
        E=double(H<0);
        properties={'square','unitary','normal','orthogonal','real', ...
            'scalable','eigenvalues'};
    otherwise
        error('butson_hadamard:Construction', ...
            'Use tao or paley.');
end
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('butson_hadamard:Mode','Invalid output mode.');
end
switch lower(char(output))
    case 'exponents'
        U=struct('E',E,'q',q,'order',n,'construction',construction);
    case {'matrix','explicit'}
        if q==2
            U=(1-2*E)/sqrt(n);
        else
            U=exp(2i*pi*E/q)/sqrt(n);
        end
    case 'spectral'
        U=unitary_spectral_reference('butson_hadamard',n);
        U.reference='Tadej and Zyczkowski (2006); Paley (1933)';
        U.notes=['Exact integer Butson phase table; normalized strict ' ...
            'unitarity. No universal analytical trace phase law. ' ...
            'Left/right Hadamard equivalence does not preserve eigenvalues.'];
    otherwise
        error('butson_hadamard:Mode','Use matrix, exponents or spectral.');
end
end
