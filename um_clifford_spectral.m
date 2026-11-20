function S=um_clifford_spectral(q,profile)
%UM_CLIFFORD_SPECTRAL Exact trace measures for elementary Clifford gates.
%   Unknown analytic laws for general gate words are explicitly distinguished
%   from the mathematical fact that their finite Clifford phases are roots
%   of unity under the fixed H/S/CNOT convention.
N=2^q;
blank=struct('kind','unknown','density',[],'support',[], ...
    'atoms',[],'moments',[]);
S=struct('family','clifford_circuit','dimension',N, ...
    'normalization','probability (density with respect to dtheta)', ...
    'finite',struct('kind','unknown','eigenvalues',[], ...
        'moments',[],'atoms',[]), ...
    'limit',blank,'ensemble',blank,'underlying',blank, ...
    'reference','Dehaene and De Moor (2003); Gross, Nezami and Walter (2021)', ...
    'notes',['The projective Clifford group is finite with the fixed ', ...
        'H/S/CNOT phase convention. No general analytic finite phase ', ...
        'multiset or CUE law is asserted for arbitrary circuit words.']);
switch profile
    case 'identity'
        z=1;w=1;
    case {'h','x','y','z'}
        z=[1;-1];w=[1/2;1/2];
    case 's'
        z=[1;1i];w=[1/2;1/2];
    case 'sdg'
        z=[1;-1i];w=[1/2;1/2];
    case {'cnot','cz','swap'}
        z=[1;-1];w=[3/4;1/4];
    otherwise
        return
end
S.finite.kind='exact';
S.finite.eigenvalues=@()repelem(z,round(N*w));
S.finite.moments=@(k)clifford_moments(z,w,k);
S.finite.atoms=struct('angles',angle(z),'weights',w);
S.notes=['Exact finite equal-weight trace eigenphase multiset for ', ...
    'a single Clifford gate, with all other qubits idle.'];
end

function m=clifford_moments(z,w,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
m=reshape(w.'*bsxfun(@power,z(:),reshape(k,1,[])),size(k));
end
