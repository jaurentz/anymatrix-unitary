function [U,properties]=daubechies_wavelet(n,levels,output)
%DAUBECHIES_WAVELET Finite periodic orthogonal four-tap (D4/db2) DWT.
%   U=DAUBECHIES_WAVELET(N,LEVELS) returns an N-by-N orthogonal
%   multilevel Daubechies four-filter-tap analysis transform, also
%   known as db2 (two vanishing moments). N must be a power of two.
%   LEVELS defaults to log2(N), and LEVELS=0 returns identity.
%   The periodized stage acts on m active averages (m even):
%      a_j = sum_(k=0)^3 h_k*x_(2j+k mod m)
%      d_j = sum_(k=0)^3 g_k*x_(2j+k mod m)
%   with h=[1+sqrt3,3+sqrt3,3-sqrt3,1-sqrt3]/(4sqrt2)
%        g=[h3,-h2,h1,-h0].
%   The stage ordering is [averages;details] and later stages act
%   only on the averages. The m=2 stage is valid by periodized
%   aliasing and reduces to a Haar-type final split.
%   Periodic boundary conventions are ESSENTIAL for unitarity.
%   OUTPUT='operator' uses O(N) two-channel filter-bank matvecs,
%   with forward/adjoint actions and multiple RHS. No claim is made
%   about generic analytic finite eigenphases of the DWT matrix.
%
%   I. Daubechies, Ten Lectures on Wavelets, SIAM, 1992.
%   https://doi.org/10.1137/1.9781611970104
properties={'square','unitary','normal','orthogonal','real', ...
            'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
maxLevels=round(log2(double(n)));
if 2^maxLevels~=n
    error('daubechies_wavelet:Order','N must be a power of two.');
end
if nargin<2||isempty(levels),levels=maxLevels;end
validateattributes(levels,{'numeric'}, ...
    {'scalar','integer','>=',0,'<=',maxLevels,'finite'});
if nargin<3||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('daubechies_wavelet:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if strcmp(output,'spectral')
    U=unitary_spectral_reference('daubechies_wavelet',n);
    U.reference='Daubechies, Ten Lectures on Wavelets (1992), D4';
    if levels==0
        U.finite.kind='exact';
        U.finite.eigenvalues=@()ones(n,1);
        U.finite.moments=@(k)ones(size(k));
        U.finite.atoms=struct('angles',0,'weights',1);
        U.notes='No decomposition levels: the transform is exactly identity.';
    else
        U.notes=['Real-orthogonal periodic D4 transform; eigenphases ' ...
            'depend on the order, level count and boundary convention.'];
    end
    return
end
if ~any(strcmp(output,{'matrix','explicit','operator','implicit'}))
    error('daubechies_wavelet:Mode','Supported modes: matrix, operator, spectral.');
end
if any(strcmp(output,{'operator','implicit'}))
    U=@(X,flag)um_daubechies_apply(X,flag,n,levels);
else
    U=um_daubechies_apply(eye(n),'notransp',n,levels);
end
end
