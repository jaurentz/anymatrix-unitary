function S=um_fibonacci_braid_spectral(n,charge,word,chirality,paths,r)
%UM_FIBONACCI_BRAID_SPECTRAL Exact finite braid-generator trace spectra.
%   The generic word has no asserted analytic finite phase multiset.
d=size(paths,1);
blank=struct('kind','unknown','density',[],'support',[], ...
    'atoms',[],'moments',[]);
S=struct('family','fibonacci_braid','dimension',d, ...
    'normalization','probability (density with respect to dtheta)', ...
    'finite',struct('kind','unknown','eigenvalues',[], ...
        'moments',[],'atoms',[]), ...
    'limit',blank,'ensemble',blank,'underlying',blank, ...
    'reference','Trebst et al. (2008); Carnahan, Zeuch and Bonesteel (2016)', ...
    'notes',['A finite Fibonacci-anyon braid representation in the ', ...
        'specified charge, fusion-tree gauge and chirality. No generic ', ...
        'analytic finite phase multiset for an arbitrary braid word.'], ...
    'anyon_count',n,'total_charge',charge,'chirality',chirality);
if isempty(word)
    z=1;counts=d;
elseif numel(word)==1
    % All adjacent braid generators are conjugate; sigma_1 is diagonal
    % in the left-associated basis with pair fusion channel a2.
    counts=[sum(paths(:,2)==0);sum(paths(:,2)==1)];
    z=r;
    if word(1)<0,z=conj(z);end
    present=counts>0;z=z(present);counts=counts(present);
else
    return
end
w=counts(:)/d;
z=z(:);
S.finite.kind='exact';
S.finite.eigenvalues=@()repelem(z,counts(:));
S.finite.atoms=struct('angles',angle(z),'weights',w);
S.finite.moments=@(k)moments(z,w,k);
S.notes=['Exact finite trace phase multiplicities for an identity ', ...
    'or a single signed braid generator, in this fusion sector.'];
end

function m=moments(z,w,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
m=reshape(w.'*bsxfun(@power,z(:),reshape(k,1,[])),size(k));
end
