function [W,properties]=electric_quantum_walk(n,omega,coinAngle,phase,boundary,mode)
%ELECTRIC_QUANTUM_WALK Finite cyclic electric two-state quantum walk.
%   W=ELECTRIC_QUANTUM_WALK(N,OMEGA,COINANGLE,PHASE,BOUNDARY,MODE)
%   constructs a 2N-by-2N strictly unitary finite-ring closure of
%   the electric quantum-walk model of Cedzich-Werner and Fan Yang.
%   Sites j=0,..,N-1 have a two-dimensional internal state in the
%   (minus,plus) order: |j,-> index 2*j+1, |j,+> index 2*j+2.
%   OMEGA is electric phase increment in CYCLES/site, PHASE is a
%   global phase in CYCLES; BOUNDARY is a twist in RADIANS. An
%   irrational real OMEGA creates a FINITE WINDOW with cyclic closure,
%   NOT a periodic irrational infinite-lattice quantum walk.
%
%   With a=cos(COINANGLE)>=0, b=sin(COINANGLE)>=0, coin
%       C=[a,b;-b,a] (order minus,plus),  C'*C=I.
%   Define the conditional shift S|j,->=|j-1,-> and
%   S|j,+>=|j+1,+>, with phase exp(-i*BOUNDARY) on the
%   minus left-wrap and exp(+i*BOUNDARY) on the plus right-wrap.
%   The electric phase D is the diagonal matrix whose two entries
%   at site j equal exp(2*pi*i*(PHASE+j*OMEGA)).
%   This generator implements EXACTLY
%       W=D*(kron(I_N,C))*S.
%
%   A second *exact* representation is W=D*L*M, where
%   L is the direct sum over sites of canonical CMV even blocks
%       Theta(b)=[b,a;a,-b],
%   and M is the cyclic odd-site swap (with the same boundary twist).
%   This is a phase-decorated CYCLIC CMV-factor model for ALL fields.
%   For OMEGA=PHASE=BOUNDARY=0, D=I and this is a genuine cyclic
%   even/odd standard-CMV product with alpha_(2j)=b,
%   alpha_(2j+1)=0. Do not identify it with the ordinary
%   finite half-line CMV generator's terminal phase.
%   For general fields, the whole-line diagonal-gauge correspondence
%   has formal skew-shift Verblunsky coefficients (Yang 2025)
%       alpha_(2j)=b*exp(2*pi*i*(j^2*OMEGA+2*j*PHASE)),
%       alpha_(2j+1)=0.
%   Those WHOLE-LINE parameters should not be confused with the
%   exact finite-ring gauge below. A finite ring always has a gauge
%   to a scalar times a twisted cyclic CMV with the SAME boundary twist;
%   absorbing the scalar is subject to a monodromy condition.
%   Explicitly, with j=0,..,N-1, g=exp(2*pi*i*(PHASE+(N-1)*OMEGA/2)):
%     alpha_(2j)=b*exp(2*pi*i*OMEGA*j*(N-1-j)), alpha_(2j+1)=0,
%     Lambda_(j,-)=exp(pi*i*OMEGA*j*(N-j)),
%     Lambda_(j,+)=exp(pi*i*OMEGA*j*(j+2-N)).
%   Then Lambda'*W*Lambda=g*L(alpha_even)*M_odd(BOUNDARY).
%   The ring electric monodromy is P=prod_j D_j=g^N.
%   For a fixed twist and positive CMV off-diagonals, an unscaled
%   representation (g=1) needs P=1 (except degenerate coins).
%
%   COINANGLE in [0,pi/2] (default .41), OMEGA default
%   (sqrt(5)-1)/2, PHASE default .13, BOUNDARY default 0.
%   N=3..100000, explicit 'matrix' and 'model' N<=2048.
%
%   MODE='matrix' (default): sparse 2N-by-2N unitary.
%   MODE='operator': sparse O(N) forward, transpose and adjoint
%       actions for arbitrary numbers of right-hand sides.
%   MODE='cmv_gauge': exact finite diagonal gauge to a GLOBAL-PHASE-
%       MULTIPLIED twisted cyclic standard CMV product. Returns sparse
%       Lambda, U_cmv, L_even, M_odd, global_phase, and Verblunsky data.
%       Lambda'*W*Lambda=U_cmv=global_phase*L_even*M_odd.
%       The cyclic boundary twist is preserved, for ANY finite OMEGA.
%       This is not an ordinary half-line finite CMV truncation.
%   MODE='factors': exact sparse L,M,D and formal skew-shift
%       Verblunsky coefficients, with finite-boundary caveat.
%   MODE='model': W,S,coin,D,L,M and original parameters.
%   MODE='spectral': exact finite phase multiset/moments when
%       COINANGLE=0 (two independent weighted N-cycles), or
%       COINANGLE=pi/2 (N independent two-state reflecting bonds).
%       Otherwise analytic finite trace law UNKNOWN; numeric eig(W)
%       is a different, non-certified reference.
%
%   At maximal reflection COINANGLE=pi/2, bond eigenvectors
%   are supported on two internal states, an exactly solvable
%   finite localization limit. This does NOT prove the irrational
%   infinite-volume Anderson localization established in [149].
%
%   [149] C. Cedzich, A. H. Werner, CMP 387 (2021), 1257-1279.
%         https://doi.org/10.1007/s00220-021-04204-w
%   [150] F. Yang, Lett. Math. Phys. 115 (2025), article 49.
%         https://doi.org/10.1007/s11005-025-01934-z
properties={'square','unitary','normal','complex', ...
    'sparse','scalable','parameter-dependent'};
if nargin==0,W=[];return;end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n) ...
        &&n==floor(n)&&n>=3&&n<=100000)
    error('electric_quantum_walk:Order','N must be an integer 3..100000.');
end
if nargin<2||isempty(omega),omega=(sqrt(5)-1)/2;end
if nargin<3||isempty(coinAngle),coinAngle=.41;end
if nargin<4||isempty(phase),phase=.13;end
if nargin<5||isempty(boundary),boundary=0;end
if nargin<6||isempty(mode),mode='matrix';end
for value={omega,coinAngle,phase,boundary}
    z=value{1};
    if ~(isnumeric(z)&&isscalar(z)&&isreal(z)&&isfinite(z))
        error('electric_quantum_walk:Parameters', ...
            'All phase, coin and boundary parameters must be real finite scalars.');
    end
end
if coinAngle<0||coinAngle>pi/2||abs(omega)>1e6||abs(phase)>1e6
    error('electric_quantum_walk:Parameters', ...
        'COINANGLE in [0,pi/2], OMEGA/PHASE bounded to 1e6 cycles.');
end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('electric_quantum_walk:Mode','MODE must be scalar text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','operator','factors','model','spectral','cmv_gauge'}))
    error('electric_quantum_walk:Mode','Unsupported output mode.');
end
if any(strcmp(mode,{'matrix','model'}))&&n>2048
    error('electric_quantum_walk:ExplicitOrder', ...
        'Explicit W and model require N<=2048; use operator/factors.');
end
j=(0:n-1)';
b=sin(coinAngle);a=cos(coinAngle);
if coinAngle==pi/2,a=0;end
Theta=sparse([b,a;a,-b]);
L=kron(speye(n),Theta);
% Odd cyclic swap: plus j -> minus j+1; minus j -> plus j-1.
colMinus=2*j+1;colPlus=2*j+2;
rowForMinus=2*mod(j-1,n)+2;
rowForPlus=2*mod(j+1,n)+1;
phaseMinus=ones(n,1);phasePlus=ones(n,1);
phaseMinus(1)=exp(-1i*boundary);
phasePlus(n)=exp(1i*boundary);
M=sparse([rowForMinus;rowForPlus], ...
         [colMinus;colPlus],[phaseMinus;phasePlus], ...
         2*n,2*n);
sitePhase=exp(2i*pi*mod(phase+omega*j,1));
d=repelem(sitePhase,2);
if strcmp(mode,'cmv_gauge')
    % Exact finite-ring identity. The scalar g accounts for the phase
    % holonomy of one circuit around the ring: g^n=prod_j sitePhase(j).
    % It cannot generally be removed without changing the CMV closure.
    g=exp(2i*pi*mod(phase+omega*(n-1)/2,1));
    monodromy=exp(2i*pi*mod(n*phase+omega*n*(n-1)/2,1));
    lambdaMinus=exp(1i*pi*mod(omega*j.*(n-j),2));
    lambdaPlus=exp(1i*pi*mod(omega*j.*(j+2-n),2));
    lambda=zeros(2*n,1);
    lambda(1:2:end)=lambdaMinus;
    lambda(2:2:end)=lambdaPlus;
    Lambda=spdiags(lambda,0,2*n,2*n);
    alphaEven=b*exp(2i*pi*mod(omega*j.*(n-1-j),1));
    Lcmv=sparse([colMinus;colMinus;colPlus;colPlus], ...
        [colMinus;colPlus;colMinus;colPlus], ...
        [conj(alphaEven);a*ones(n,1);a*ones(n,1);-alphaEven], ...
        2*n,2*n);
    Ucmv=g*Lcmv*M;
    W=struct('Lambda',Lambda,'lambda_diagonal',lambda, ...
        'global_phase',g,'monodromy',monodromy, ...
        'alpha_even',alphaEven,'alpha_odd',zeros(n,1), ...
        'L_even',Lcmv,'M_odd',M,'U_cmv',Ucmv, ...
        'boundary_twist',boundary, ...
        'notes',['Exact finite-ring gauge: Lambda''*W*Lambda = ', ...
          'U_cmv = global_phase*L_even*M_odd for every finite ', ...
          'field and twist; this is cyclic, not half-line CMV. ', ...
          'The formal whole-line skew-shift coefficients in ', ...
          'factors mode are a different gauge convention. ', ...
          'Unscaled fixed-twist CMV needs monodromy 1 ', ...
          'for nondegenerate coin.']);
    return
end
if strcmp(mode,'operator')
    W=@(X,flag)um_electric_walk_apply(X,flag,d,L,M);
    return
end
if strcmp(mode,'spectral')
    blank=struct('kind','unknown','density',[], ...
       'support',[],'atoms',[],'moments',[]);
    finite=struct('kind','unknown','eigenvalues',[], ...
       'moments',[],'atoms',[]);
    description=['Generic finite electric quantum walk eigenphases are ', ...
      'not assigned a closed analytic trace law. Finite cyclic ', ...
      'irrational-field closure is NOT the infinite irrational ', ...
      'lattice; localization and full-circle spectrum theorems ', ...
      'cannot be inferred from one finite W.'];
    if coinAngle==0
        total=2*pi*(n*phase+omega*n*(n-1)/2);
        roots=2*pi*(0:n-1)'/n;
        z=[exp(1i*((total-boundary)/n+roots)); ...
           exp(1i*((total+boundary)/n+roots))];
        finite.kind='exact';
        finite.eigenvalues=@()z;
        finite.moments=@(k)um_electric_walk_moments(z,k);
        description=['Exact transparent-coin limit: independent ', ...
           'minus and plus N-site electric weighted cycles with ', ...
           'opposite boundary twist.'];
    elseif coinAngle==pi/2
        next=mod(j+1,n)+1;
        center=(angle(sitePhase)+angle(sitePhase(next)))/2;
        z=[1i*exp(1i*center);-1i*exp(1i*center)];
        finite.kind='exact';
        finite.eigenvalues=@()z;
        finite.moments=@(k)um_electric_walk_moments(z,k);
        description=['Exact fully reflecting-coin limit: N decoupled ', ...
           'two-state bond eigenspaces, each with phase pair ', ...
           '+/-i exp(i*(phi_j+phi_(j+1))/2); boundary twist cancels ', ...
           'within a bond. This is finite bond localization only.'];
    end
    W=struct('family','electric_quantum_walk', ...
      'dimension',2*n, ...
      'normalization','probability (density with respect to dtheta)', ...
      'finite',finite,'limit',blank,'ensemble',blank, ...
      'underlying',blank, ...
      'reference','Cedzich-Werner (2021); Yang (2025)', ...
      'notes',description,'field',omega,'coin_angle',coinAngle, ...
      'boundary_twist',boundary);
    return
end
D=spdiags(d,0,2*n,2*n);
if strcmp(mode,'factors')
    formalEven=b*exp(2i*pi*mod(j.^2*omega+2*j*phase,1));
    W=struct('D',D,'L_even',L,'M_odd',M, ...
       'field',omega,'phase',phase,'boundary_twist',boundary, ...
       'coin',full([a,b;-b,a]),'formal_alpha_even',formalEven, ...
       'formal_alpha_odd',zeros(n,1), ...
       'standard_cyclic_cmv',omega==0&&phase==0&&boundary==0, ...
       'notes',['W = D*L_even*M_odd exactly on this ring. ', ...
         'The formal skew-shift coefficients describe the ', ...
         'infinite-lattice diagonal gauge, not necessarily ', ...
         'the same periodic finite boundary condition.']);
    return
end
W=D*L*M;
if strcmp(mode,'model')
    F=kron(speye(n),sparse([0 1;1 0]));
    shift=F*M;
    C=kron(speye(n),sparse([a,b;-b,a]));
    W=struct('U',W,'D',D,'S',shift,'coin',C, ...
      'L_even',L,'M_odd',M,'site_phase',sitePhase, ...
      'coin_angle',coinAngle,'field',omega, ...
      'phase',phase,'boundary_twist',boundary, ...
      'notes','Exact finite W=D*coin*S=D*L_even*M_odd; S=F*M_odd.');
end
end

function moments=um_electric_walk_moments(z,k)
validateattributes(k,{'numeric'},{'real','finite','integer'});
moments=reshape(mean(bsxfun(@power,z(:),reshape(k,1,[])),1),size(k));
end
