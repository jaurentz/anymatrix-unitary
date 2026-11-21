function [U,properties]=fourier_zernike_disk(nRad,nPhi,mode)
%FOURIER_ZERNIKE_DISK Exact finite unitary modal Zernike/Fourier hybrid.
%   The positive-leading radial polynomials are the classical
%     R_{m+2k}^{|m|}(r)=r^|m| P_k^(0,|m|)(2r^2-1).
%   This generator combines nRad-by-nRad orthogonal Gaussian
%   radial transforms Q_m with an orthonormal nPhi-point angular
%   INVERSE DFT, acting on an nRad-by-nPhi coefficient array C(k,m).
%
%   INPUT is vectorized with RADIAL INDEX FASTEST, Fourier mode index
%   slowest (MATLAB column-major). FFT columns k=0,...,nPhi-1
%   correspond to m=k for k<=floor(nPhi/2), else m=k-nPhi.
%   For each mode m, Gaussian nodes r_(j,m) are the square roots
%   of (x_j+1)/2 with x_j nRad-point Gauss-Jacobi alpha=0,beta=|m|.
%   The radial matrix Q_m(j,k) is the normalized quadrature-weighted
%   evaluation of P_(k-1)^(0,|m|)(x_j) and is exactly real orthogonal
%   in mathematical arithmetic, by Gauss exactness through degree
%   2*nRad-1. U = (F_phi^* tensor I_rad)*blockdiag(Q_m),
%   with F_phi normalized DFT. Thus U'*U=I EXACTLY in real arithmetic.
%
%   The output angular synthesis at a given radial SLOT j involves
%   DIFFERENT PHYSICAL RADII r_(j,m) for different modes m. Hence
%   this output is a MODAL-RADIAL-HYBRID ARRAY, NOT samples on a
%   common tensor disk grid. No square common-grid collocation
%   unitarity is asserted. Standard physical disk interpolation
%   would require additional NONUNITARY radial interpolation.
%
%   MODE='matrix' (default): dense U for N=nRad*nPhi<=512.
%       'operator': precomputed radial Q_m plus FFT with
%                   O(nPhi*nRad^2+nRad*nPhi*log(nPhi)) work/RHS,
%                   O(nPhi*nRad^2) block storage after O(nPhi*nRad^3)
%                   dense Golub-Welsch setup. Forward/adjoint/
%                   ordinary transpose, complex multi-RHS, size.
%       'model': radial nodes and normalized mode-specific Gaussian
%                weights, radial matrices, signed m modes, degree
%                labels, explicit warning about non-common grid.
%       'spectral': generic finite matrix eigenphases UNKNOWN;
%                   radial quadrature nodes are NOT eigenphases.
%
%   The transform is nontrivial and application-derived: for each
%   |m|, orthonormal disk radial Zernike functions are
%      sqrt(2*(2*k+|m|+1))*r^|m| P_k^(0,|m|)(2r^2-1).
%   Their norm is taken with r dr and the normalized angular
%   Fourier factors. A mode-dependent Gaussian rule for
%   r^(2*|m|+1) dr makes the finite radial block orthogonal.
%   The physical radial quadrature weights wPhysical are
%   wJacobiNormalized/[2*(|m|+1)*r^(2*|m|)] and therefore vary
%   between Fourier modes: NEVER interpret them as a common
%   quadrature rule for arbitrary disk functions.
%
%   nRad<=128, 2<=nPhi<=64, nRad*nPhi<=4096 for operator;
%   no genuinely fast arbitrary Jacobi radial transform is claimed.
%   [174] Janssen (2014) DOI 10.1364/JOSAA.31.001604.
%   [175] Slevinsky, FastTransforms Zernike/disk, official documentation.
%         https://mikaelslevinsky.github.io/FastTransforms/transforms.html
properties={'square','unitary','normal','complex','scalable', ...
    'parameter-dependent'};
if nargin==0,U=[];return;end
if ~(isnumeric(nRad)&&isscalar(nRad)&&isreal(nRad) ...
        &&isfinite(nRad)&&nRad==fix(nRad)&&nRad>=1&&nRad<=128)
    error('fourier_zernike_disk:RadialOrder','Require 1<=nRad<=128 integer.');
end
if nargin<2||isempty(nPhi),nPhi=7;end
if ~(isnumeric(nPhi)&&isscalar(nPhi)&&isreal(nPhi) ...
        &&isfinite(nPhi)&&nPhi==fix(nPhi)&&nPhi>=2&&nPhi<=64)
    error('fourier_zernike_disk:AngularOrder','Require 2<=nPhi<=64 integer.');
end
N=nRad*nPhi;
if N>4096
    error('fourier_zernike_disk:Order','Operator order must be <=4096.');
end
if nargin<3||isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('fourier_zernike_disk:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','operator','model','spectral'}))
    error('fourier_zernike_disk:Mode', ...
        'Use matrix, operator, model, spectral.');
end
k=0:nPhi-1;
m=k;
m(k>floor(nPhi/2))=m(k>floor(nPhi/2))-nPhi;
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('fourier_zernike_disk',N);
    U.radial_order=nRad;U.angular_order=nPhi;
    U.signed_modes=m;
    U.reference='Janssen 2014; Slevinsky FastTransforms Zernike';
    U.notes=['Finite matrix eigenphases unknown. Gauss-Jacobi ', ...
        'nodes depend on azimuthal order and are NOT unit-circle ', ...
        'eigenphases; output is modal-radius hybrid, NOT common ', ...
        'physical grid samples.'];
    return
end
if strcmp(mode,'matrix')&&N>512
    error('fourier_zernike_disk:DenseOrder', ...
        'Dense matrix output requires nRad*nPhi<=512.');
end
Q=cell(nPhi,1);
radialNodes=zeros(nRad,nPhi);
radialWeights=zeros(nRad,nPhi);
degree=zeros(nRad,nPhi);
distinct=unique(abs(m),'stable');
for t=1:numel(distinct)
    a=distinct(t);
    diagonal=zeros(nRad,1);
    diagonal(1)=a/(a+2);
    if nRad>1
        j=(1:nRad-1)';
        diagonal(2:end)=a^2./((2*j+a).*(2*j+a+2));
        off=zeros(nRad-1,1);
        off(1)=2/(a+2)*sqrt((a+1)/(a+3));
        if nRad>2
            jj=(2:nRad-1)';
            off(2:end)=2*jj.*(jj+a)./ ...
                ((2*jj+a).*sqrt((2*jj+a-1).*(2*jj+a+1)));
        end
    else
        off=[];
    end
    [Qone,info]=um_gauss_basis(diagonal,off);
    r=sqrt((info.nodes+1)/2);
    % For normalized Gauss-Jacobi weight (1+x)^a, the physical
    % radial discrete measure for r dr includes the factor r^(-2a).
    physicalWeights=info.normalized_weights./ ...
        (2*(a+1)*r.^(2*a));
    columns=find(abs(m)==a);
    for col=columns
        Q{col}=Qone;
        radialNodes(:,col)=r;
        radialWeights(:,col)=physicalWeights;
        degree(:,col)=(a+2*(0:nRad-1))';
    end
end
apply=@(X,flag)um_fourier_zernike_apply(X,flag,Q,nRad,nPhi);
if strcmp(mode,'operator'),U=apply;return;end
if strcmp(mode,'matrix')
    U=apply(eye(N),'notransp');
    return
end
U=struct('dimension',N,'radial_order',nRad, ...
    'angular_order',nPhi,'signed_modes',m(:), ...
    'radial_nodes_by_mode',radialNodes, ...
    'physical_radial_weights_by_mode',radialWeights, ...
    'radial_degree_by_mode',degree, ...
    'radial_blocks',{Q},'angular_angles',2*pi*(0:nPhi-1)/nPhi, ...
    'apply',apply,'is_common_physical_grid',false, ...
    'notes',['Each azimuthal mode m has its own Gauss-Jacobi ', ...
    'radial nodes. The U output is an angularly synthesized ', ...
    'modal-radial array, NOT physical disk samples at a shared ', ...
    'set of radii; physical radial weights are mode-specific. ', ...
    'Mathematically exact finite unitarity is by direct sum of ', ...
    'positive Gauss quadrature blocks and a unitary Fourier map.']);
end
