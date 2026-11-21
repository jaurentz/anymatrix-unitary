function [U,properties]=so3_peter_weyl_transform(L,m,n,output)
%SO3_PETER_WEYL_TRANSFORM Strict square fixed-(m,n) Wigner-Jacobi block.
%   U=SO3_PETER_WEYL_TRANSFORM(L,M,N) is the REAL ORTHOGONAL square
%   P-by-P fixed pair of Wigner D Euler-index modes, with
%      P=L-max(abs(M),abs(N)) >0, ell=max(abs(M),abs(N)),...,L-1.
%   The normalized Haar measure on SO(3) in z-y-z Euler angles is
%      dmu=dalpha * sin(beta)*dbeta * dgamma/(8*pi^2),
%      alpha,gamma in [0,2*pi), beta in [0,pi].
%   Full Peter--Weyl basis functions are
%      F_(ell,m,n)(alpha,beta,gamma)=sqrt(2ell+1)
%                  *exp(-i*m*alpha)*d^ell_(m,n)(beta)*exp(-i*n*gamma).
%
%   For a fixed pair (m,n), set j0=max(abs(m),abs(n)),
%   a=abs(m-n), b=abs(m+n), p=ell-j0.
%   The latitude basis G_ell(x)=sqrt((2ell+1)/2)
%   *d^ell_(m,n)(acos(x)) equals a constant signed Jacobi
%   polynomial times (1-x)^(a/2)*(1+x)^(b/2).
%   The square normalized P-node Gauss-Jacobi evaluation matrix
%   U=(-1)^max(0,m-n)*GAUSS_JACOBI_TRANSFORM(P,a,b) is orthogonal
%   exactly in mathematical arithmetic. Its eigenphases are
%   UNKNOWN; the real quadrature nodes are NOT its eigenvalues.
%
%   OUTPUT='matrix' (default): strict square real orthogonal block.
%          'spectral': unknown finite trace phase law.
%          'parameters': dimensions, normalized Haar and quadrature.
%          'sampling': RECTANGULAR full bandlimited SO(3) weighted
%                 evaluation A, with coefficients sum_(ell=0..L-1)
%                 (2ell+1)^2=L*(4L^2-1)/3, and sampling rows
%                 L*(2L-1)^2. The nodes are L-point Gauss-Legendre
%                 beta cosine x_j and Mphi=2L-1 equispaced
%                 alpha/gamma angles, with normalized HAAR weights
%                 wbeta_j/(2*Mphi^2). Each column is
%                    sqrt(wbeta_j/(2*Mphi^2)) *
%                     sqrt(2ell+1)*D^ell_(m,n)(alpha,beta,gamma).
%                 Consequently A'*A=I (exact in arithmetic),
%                 although A is NOT a square unitary unless L=1.
%                 Row order beta fastest, then alpha, then gamma;
%                 column order ell increasing, m=-ell:ell,
%                 n=-ell:ell.
%          'sampling_operator': exact FFT-separated rectangular
%                 synthesis/adjoint using an O(L^4)-storage cache
%                 of Wigner-Jacobi latitude blocks and two FFTs.
%                 Per RHS O(L^4+L^3 log L) arithmetic rather than
%                 O(L^6) dense full SO(3) evaluation.
%
%   Fixed blocks and sampling are different from WIGNER_ROTATION:
%   that family returns the irreducible D^ell(R) which acts
%   on a single coefficient vector, while this generator samples
%   bandlimited scalar functions on the entire group SO(3).
%
%   Explicit rectangular 'sampling' capped at L<=7;
%   'sampling_operator' and square block support L<=20.
%   FFT separation is an independent Gauss-Legendre realization,
%   NOT a claim to implement the original MW torus-extension
%   sampling method or every optimized S2FFT/Wigner recurrence.
%
%   [166] McEwen et al. A Novel Sampling Theorem on SO(3), 2015:
%   https://doi.org/10.1109/LSP.2015.2490676
%   [167] Khalid et al. Gauss-Legendre Sampling on SO(3), 2016:
%   https://doi.org/10.1109/LSP.2015.2503295
%   [108] Gimbutas-Greengard Wigner rotations (2009):
%   https://doi.org/10.1016/j.jcp.2009.05.014
%   [104] Golub-Welsch (1969):
%   https://doi.org/10.1090/S0025-5718-69-99647-1
properties={'square','unitary','orthogonal','real','normal', ...
    'scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'}, ...
    {'scalar','integer','positive','finite','<=',20});
if nargin<2||isempty(m),m=0;end
if nargin<3||isempty(n),n=0;end
validateattributes(m,{'numeric'},{'scalar','integer','real','finite'});
validateattributes(n,{'numeric'},{'scalar','integer','real','finite'});
if max(abs(m),abs(n))>=L
    error('so3_peter_weyl_transform:Order', ...
        'Require max(abs(m),abs(n))<L.');
end
if nargin<4||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('so3_peter_weyl_transform:Mode','OUTPUT must be scalar text.');
end
output=lower(char(output));
if ~any(strcmp(output, ...
        {'matrix','spectral','parameters','sampling','sampling_operator'}))
    error('so3_peter_weyl_transform:Mode','Unknown output mode.');
end
j0=max(abs(m),abs(n));P=L-j0;
a=abs(m-n);b=abs(m+n);phase=(-1)^max(0,m-n);
if strcmp(output,'matrix')
    U=phase*gauss_jacobi_transform(P,a,b);
    return
end
if strcmp(output,'spectral')
    U=unitary_spectral_reference('so3_peter_weyl_transform',P);
    U.bandlimit=L;U.left_frequency=m;U.right_frequency=n;
    U.reference='Khalid et al., IEEE SPL 23 (2016), 207-211';
    U.notes=['Exact square fixed-Wigner-index Gauss-Jacobi latitude ', ...
        'evaluation with unknown finite matrix eigenphases. The ', ...
        'full SO(3) Peter-Weyl transform is rectangular isometric, ', ...
        'NOT a square unitary or an irreducible rotation matrix.'];
    return
end
[x,w]=um_legendre_gauss_rule(L);M=2*L-1;
C=L*(4*L^2-1)/3;rows=L*M^2;
if strcmp(output,'parameters')
    U=struct('bandlimit',L,'fixed_m',m,'fixed_n',n, ...
        'fixed_block_dimension',P,'fixed_jacobi_alpha',a, ...
        'fixed_jacobi_beta',b,'fixed_block_sign',phase, ...
        'beta_cos_nodes',x,'beta_weights_dx',w, ...
        'alpha_nodes',2*pi*(0:M-1)/M, ...
        'gamma_nodes',2*pi*(0:M-1)/M, ...
        'normalized_haar_weights',w/(2*M^2), ...
        'sampling_rows',rows,'sampling_columns',C, ...
        'row_order','beta fastest, then alpha, then gamma', ...
        'coefficient_order','ell, m=-ell:ell, n=-ell:ell', ...
        'convention','sqrt(2ell+1)*exp(-im alpha)*d_(m,n)(beta)*exp(-in gamma)');
    return
end
if strcmp(output,'sampling')&&L>7
    error('so3_peter_weyl_transform:SamplingSize', ...
        'Explicit rectangular sampling limited to L<=7; use sampling_operator.');
end
data=struct('L',L,'M',M,'Ncoef',C);
data.blocks=cell(M,M);data.cols=cell(M,M);
for mm=-(L-1):L-1
    i=mod(mm,M)+1;
    for nn=-(L-1):L-1
        j=mod(nn,M)+1;
        G=um_wigner_jacobi(L,mm,nn,x);
        data.blocks{i,j}=bsxfun(@times,sqrt(w),G);
    end
end
index=0;
for ell=0:L-1
    for mm=-ell:ell
        i=mod(mm,M)+1;
        for nn=-ell:ell
            index=index+1;j=mod(nn,M)+1;
            data.cols{i,j}=[data.cols{i,j} index]; %#ok<AGROW>
        end
    end
end
if strcmp(output,'sampling_operator')
    U=@(X,flag)um_so3_sampling_apply(X,flag,data);
else
    F=@(X,flag)um_so3_sampling_apply(X,flag,data);
    U=F(eye(C),'notransp');
end
end
