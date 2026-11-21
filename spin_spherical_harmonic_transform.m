function [U,properties]=spin_spherical_harmonic_transform(L,s,m,output)
%SPIN_SPHERICAL_HARMONIC_TRANSFORM Square fixed-order spin-s latitude Q.
%   U=SPIN_SPHERICAL_HARMONIC_TRANSFORM(L,S,M) constructs a genuine
%   REAL ORTHOGONAL (L-j0)-square fixed azimuth m / spin s block,
%   j0=max(abs(s),abs(m)), ell=j0,...,L-1, L>j0.
%   The Wigner d-element d^ell_(m,-s)(theta) equals
%      (-1)^max(0,m+s) * sqrt(2/(2ell+1)) *
%      (1-x)^(a/2)*(1+x)^(b/2) * p_(ell-j0)(x)/sqrt(mu0),
%   x=cos(theta), a=abs(m+s), b=abs(m-s), with normalized,
%   positive-leading Jacobi p and mu0=int_(−1)^1 (1-x)^a(1+x)^b dx.
%   We use the spin harmonics
%      {}_sY_(ell,m)(theta,phi)
%         =(-1)^s*sqrt((2ell+1)/(4pi))*d^ell_(m,-s)(theta)*exp(i*m*phi).
%   At the (L-j0)-node Gauss--Jacobi rule (a,b), the exact weighted
%   latitude block is U=(-1)^(s+max(0,m+s))*GAUSS_JACOBI_TRANSFORM(
%   L-j0,a,b). It is square and orthogonal in mathematical arithmetic;
%   its eigenphases are UNKNOWN (the real Jacobi zeros are not them).
%   Default S=0,M=0 recovers the ordinary Gauss-Legendre block.
%
%   OUTPUT='matrix' (default) : square real orthogonal fixed (s,m).
%          'spectral' : finite transform phase law UNKNOWN.
%          'parameters': exact spin/Jacobi conventions; complete
%                        full-sphere sampling grid and coefficient labels.
%          'sampling': full-sphere RECTANGULAR weighted evaluation A,
%                      size [L*(2L-1),L^2-s^2], with L Gauss-Legendre
%                      x=cos(theta) nodes and Mphi=2L-1 equally spaced
%                      azimuths 2pi*k/Mphi. Rows: theta FASTEST, phi
%                      SLOWEST; columns: ell increasing, m=-ell:ell.
%                      Sample weights wtheta_j*(2pi/Mphi).
%                      A'*A=I EXACTLY by degree<=2L-2 polynomial
%                      quadrature and exact finite azimuth Fourier
%                      orthogonality. It is NOT a square unitary unless
%                      L=1,S=0.
%          'sampling_operator': rectangular FFT-separated A and A'
%                      with O(L^3) precomputation/storage for Wigner
%                      latitude blocks and O(L^3+L^2 log L) apply,
%                      vs O(L^4) dense A multiplication. This is an
%                      exact bandlimited GL sampling formulation, NOT
%                      the McEwen-Wiaux equiangular sampling grid.
%
%   Uses integer spin |S|<L, integer azimuth |M|<L; permits S=+-1,
%   S=+-2 and the associated conjugation relation
%      conj({}_sY_(ell,m))=(-1)^(m+s)*{}_-sY_(ell,-m).
%   Spin sample frames have a polar coordinate singularity; no
%   pointwise pole values or global smooth spin frame are claimed.
%   Dense 'sampling' capped L<=12, fast FFT sampling operator L<=32.
%   The default square block reuses an O(N^3) dense Gauss-Jacobi
%   eigensolver; it is NOT FFT-fast at generic (a,b).
%
%   [165] Huffenberger-Wandelt, ApJS 189 (2010), 255-260:
%   https://doi.org/10.1088/0067-0049/189/2/255
%   [107] Slevinsky, Fast backward stable Fourier conversions:
%   https://arxiv.org/abs/1705.05448
%   [104] Golub-Welsch: https://doi.org/10.1090/S0025-5718-69-99647-1
properties={'square','unitary','normal','orthogonal','real', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',32});
if nargin<2||isempty(s),s=0;end
if nargin<3||isempty(m),m=0;end
validateattributes(s,{'numeric'},{'scalar','integer','real','finite'});
validateattributes(m,{'numeric'},{'scalar','integer','real','finite'});
if max(abs(s),abs(m))>=L
    error('spin_spherical_harmonic_transform:Order', ...
        'Need max(abs(spin),abs(azimuth))<L.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('spin_spherical_harmonic_transform:Mode', ...
        'OUTPUT must be scalar text.');
end
output=lower(char(output));
if ~any(strcmp(output, ...
        {'matrix','spectral','parameters','sampling','sampling_operator'}))
    error('spin_spherical_harmonic_transform:Mode','Unknown output mode.');
end
j0=max(abs(s),abs(m));N=L-j0;
a=abs(m+s);b=abs(m-s);
signPhase=(-1)^(s+max(0,m+s));
if strcmp(output,'matrix')
    U=signPhase*gauss_jacobi_transform(N,a,b);
    return
end
if strcmp(output,'spectral')
    U=unitary_spectral_reference('spin_spherical_harmonic_transform',N);
    U.spin=s;U.azimuth=m;U.bandlimit=L;
    U.reference='Huffenberger-Wandelt, ApJS 189 (2010), 255-260';
    U.notes=['Fixed-(spin,m) square real Gauss-Jacobi latitude transform. ', ...
        'The finite phase spectrum is unknown; real Gauss zeros are ', ...
        'quadrature coordinates, not eigenphases. The full-sphere ', ...
        'operator is a separate rectangular isometry.'];
    return
end
[x,w]=um_legendre_gauss_rule(L);
M=2*L-1;
if strcmp(output,'parameters')
    U=struct('bandlimit',L,'spin',s,'azimuth',m,'ell_min',j0, ...
        'fixed_block_dimension',N,'jacobi_alpha',a,'jacobi_beta',b, ...
        'fixed_block_sign',signPhase,'theta_cos_nodes',x, ...
        'theta_weights',w,'phi_nodes',2*pi*(0:M-1)/M, ...
        'sphere_sample_weights',w*(2*pi/M), ...
        'sampling_rows',L*M,'sampling_columns',L^2-s^2, ...
        'coefficient_order','ell increasing; m=-ell:ell', ...
        'sampling_row_order','theta fastest; phi slowest', ...
        'spin_harmonic_convention', ...
        '(-1)^s sqrt((2ell+1)/(4pi)) d_(m,-s) exp(i*m*phi)');
    return
end
if strcmp(output,'sampling')&&L>12
    error('spin_spherical_harmonic_transform:SamplingSize', ...
        'Explicit rectangular sampling capped L<=12; use sampling_operator.');
end
data=struct('L',L,'M',M,'Ncoef',L^2-s^2);
data.blocks=cell(2*L-1,1);
data.cols=cell(2*L-1,1);
% The phase (-1)^s and sign of d are retained in every block.
for az=-(L-1):(L-1)
    j=az+L;
    if max(abs(az),abs(s))>=L,continue;end
    G=um_wigner_jacobi(L,az,-s,x);
    data.blocks{j}=(-1)^s*bsxfun(@times,sqrt(w),G);
end
labels=zeros(L^2-s^2,2);col=0;
for ell=abs(s):L-1
    for az=-ell:ell
        col=col+1;labels(col,:)=[ell az];
        data.cols{az+L}=[data.cols{az+L} col]; %#ok<AGROW>
    end
end
if strcmp(output,'sampling_operator')
    U=@(X,flag)um_spin_sampling_apply(X,flag,data);
else
    F=@(X,flag)um_spin_sampling_apply(X,flag,data);
    U=F(eye(data.Ncoef),'notransp');
end
end
