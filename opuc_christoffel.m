function [U,info]=opuc_christoffel(n,nodes,weights,zeta,representation,output)
%OPUC_CHRISTOFFEL Positive Christoffel change of a finite circle measure.
%   [U,INFO]=OPUC_CHRISTOFFEL(N,NODES,WEIGHTS,ZETA,REP,OUTPUT)
%   applies dmu_new(z)=|z-ZETA|^2*dmu(z)/Z to the strictly
%   positive finite atomic measure dmu=sum_j w_j delta(z_j),
%   with Z=sum_j w_j*|z_j-ZETA|^2 > 0.
%   A support atom is removed when ZETA equals that node exactly;
%   then the resulting unitary has reduced order. All remaining weights
%   are strictly positive, so the finite inverse Szego construction
%   (INVERSE_OPUC) gives a unitary CMV or Hessenberg model.
%   When all support atoms remain, the new matrix has the SAME
%   finite trace eigenphase multiset but DIFFERENT e1 cyclic weights.
%   This is a Christoffel change of MEASURE, not a CMV Darboux
%   matrix factorization, and no general isospectrality is claimed
%   for support deletion or infinite-dimensional measures.
%   REP='cmv' (default), 'hessenberg', or 'inverse_hessenberg'.
%   OUTPUT='matrix' (default), 'operator', 'parameters', 'spectral'.
%   INFO provides original and transformed support, normalized weights,
%   Christoffel normalization Z, removed atom indices, and dimensions.
%
%   [5] B. Simon, OPUC: https://doi.org/10.1090/coll054.1
%   [Darboux] Cantero et al.: https://doi.org/10.1016/j.aim.2016.03.042
if nargin<1,error('opuc_christoffel:Arguments', ...
        'Supply N, circle nodes and weights.');end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(nodes)
    nodes=exp(2i*pi*((0:n-1)'+.17)/n);
end
if nargin<3||isempty(weights),weights=ones(n,1)/n;end
if nargin<4||isempty(zeta),zeta=.3+.2i;end
if ~isnumeric(zeta)||~isscalar(zeta)||~isfinite(zeta)
    error('opuc_christoffel:Zeta','ZETA must be finite complex scalar.');
end
if nargin<5||isempty(representation),representation='cmv';end
if nargin<6||isempty(output),output='matrix';end
% Delegate strict support/weight validation to the existing inverse
% spectral routine, with no duplicated Verblunsky implementation.
old=inverse_opuc(n,nodes,weights,'cmv','parameters');
z=old.nodes;
w=old.weights;
factor=abs(z-zeta).^2;
if any(~isfinite(factor))
    error('opuc_christoffel:Zeta','Christoffel factors overflowed.');
end
normalization=sum(w.*factor);
if ~isfinite(normalization)||normalization<=0
    error('opuc_christoffel:ZeroMeasure', ...
        'The Christoffel factor annihilates the entire measure.');
end
keep=factor>0;
newNodes=z(keep);
newWeights=w(keep).*factor(keep)/normalization;
% A factor can underflow or generate extremely small weights; the
% underlying inverse recurrence rejects near-breakdown cases.
if any(~isfinite(newWeights))||any(newWeights<=0)
    error('opuc_christoffel:WeightUnderflow', ...
        'Transformed positive weights cannot be represented.');
end
newWeights=newWeights/sum(newWeights);
info=struct('original_dimension',n,'dimension',numel(newNodes), ...
    'original_nodes',z,'original_weights',w, ...
    'nodes',newNodes,'weights',newWeights, ...
    'removed_indices',find(~keep),'normalization',normalization, ...
    'factors',factor,'zeta',zeta, ...
    'transformation','finite_positive_Christoffel_measure_change');
U=inverse_opuc(numel(newNodes),newNodes,newWeights,representation,output);
if (ischar(output)||isstring(output))&&strcmpi(char(output),'spectral')
    U.christoffel=info;
    U.notes=[U.notes ' Christoffel reweighting preserves the finite ' ...
        'trace phases only while all original support nodes remain.'];
end
end
