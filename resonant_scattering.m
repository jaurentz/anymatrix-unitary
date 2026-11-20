function [S,properties]=resonant_scattering(n,E,H,W,mode)
%RESONANT_SCATTERING Lossless flux-normalized multichannel scattering.
%   S=RESONANT_SCATTERING(N,E,H,W,MODE) is the N-by-N unitary
%   on-shell scattering matrix for a wide-band, finite-dimensional,
%   Hermitian internal Hamiltonian H (M-by-M) coupled through W
%   (M-by-N) to N orthonormal, flux-normalized open asymptotic channels.
%   All open channels are included; no absorption, channel thresholds,
%   or energy-dependent coupling. There are M<=N internal states.
%
%   Gamma=W*W', A(z)=z*I-H+i*Gamma, and
%        S(E)=I-2*i*W'*(A(E)\W),  E real.
%   H must be Hermitian and W full ROW rank (Gamma positive definite).
%   This guarantees A(E) invertible for every real E and resonance
%   poles eig(H-i*Gamma) strictly below the real energy axis.
%   Pole locations are NOT the internal Hermitian eigenvalues eig(H).
%   For real H,W the scattering matrix is complex symmetric (TRS).
%   Default H is a reproducible real-symmetric tridiagonal cavity,
%   default W diagonal with positive channel coupling rates.
%
%   MODE:
%     'matrix' (default): return N-by-N scattering S(E).
%     'time_delay': Hermitian positive semidefinite Wigner-Smith
%       matrix Q(E)=-i*S(E)'*dS/dE=2*(A(E)\W)'*(A(E)\W)
%       in units of hbar=1. This is NOT a unitary generator output.
%     'derivative': analytic dS(E)/dE, NOT unitary.
%     'model': return H, W, Gamma, effective Hamiltonian, resonance
%       poles, and function handles S(z), derivative and time_delay.
%       The meromorphic S(z) is unitary ONLY for real z.
%     'spectral': return finite trace-measure reference data. Exact
%       finite phases are available for diagonal H and orthogonal
%       coupling rows (diagonal W*W'); other cases are marked UNKNOWN.
%   Do not equate those finite phase multiset(s), Wigner-Smith
%   proper delay times, and the complex resonance poles.
%
%   References:
%   [132] C. W. J. Beenakker, Rev. Mod. Phys. 69 (1997), 731-808.
%   https://doi.org/10.1103/RevModPhys.69.731
%   [133] P. W. Brouwer, K. M. Frahm, C. W. J. Beenakker,
%   Phys. Rev. Lett. 78 (1997), 4737-4740.
%   https://doi.org/10.1103/PhysRevLett.78.4737
%   [134] L. Erdos, T. Kruger, Y. Nemish, Ann. Henri Poincare
%   22 (2021), 4205-4269.
%   https://doi.org/10.1007/s00023-021-01085-6
properties={'square','unitary','normal','complex','parameter-dependent'};
if nargin==0,S=[];return;end
if ~(isnumeric(n)&&isscalar(n)&&isreal(n)&&isfinite(n) ...
        &&n==floor(n)&&n>=1&&n<=512)
    error('resonant_scattering:Channels','N must be an integer from 1 to 512.');
end
if nargin<2 || isempty(E),E=0.2;end
if nargin<3 || isempty(H)
    d=(0:n-1)';
    H=diag(.65*(d-(n-1)/2));
    if n>1
        off=.21+.06*sin((1:n-1)');
        H=H+diag(off,1)+diag(off,-1);
    end
end
if nargin<4 || isempty(W)
    m=size(H,1);
    if m>n
        error('resonant_scattering:Channels', ...
            'Number of open channels must be at least internal order.');
    end
    W=diag(sqrt(.35+.09*(1:m)'/m));
    if m<n,W=[W,zeros(m,n-m)];end
end
if nargin<5 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('resonant_scattering:Mode','Mode must be text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','time_delay','derivative','model','spectral'}))
    error('resonant_scattering:Mode','Invalid output mode.');
end
validateattributes(H,{'numeric'},{'2d','square','finite','nonempty'});
m=size(H,1);
if m>n
    error('resonant_scattering:Channels', ...
        'At least M open channels are required for M internal states.');
end
if ~(isnumeric(W)&&ismatrix(W)&&isequal(size(W),[m n]) ...
        &&all(isfinite(W(:))))
    error('resonant_scattering:CouplingShape', ...
        'W must be a finite M-by-N coupling matrix.');
end
if norm(H-H','fro')>1e-11*max(1,norm(H,'fro'))
    error('resonant_scattering:Hamiltonian','H must be Hermitian.');
end
H=(H+H')/2;
Gamma=W*W';
[~,p]=chol((Gamma+Gamma')/2);
if p~=0 || rcond(Gamma)<1e-12
    error('resonant_scattering:Coupling', ...
        'W must have full row rank with well-conditioned W*W''.');
end
if strcmp(mode,'model')
    heff=H-1i*Gamma;
    S=struct('dimension',n,'internal_order',m, ...
        'H',H,'W',W,'Gamma',Gamma, ...
        'effective_hamiltonian',heff, ...
        'poles',eig(heff), ...
        'S',@(z)um_resonant_scattering_eval(z,H,W,'matrix'), ...
        'derivative',@(z)um_resonant_scattering_eval(z,H,W,'derivative'), ...
        'time_delay',@(z)um_resonant_scattering_eval(z,H,W,'time_delay'), ...
        'channel_normalization','orthonormal incoming/outgoing unit-flux channels', ...
        'notes',['S(z) is meromorphic, on-shell unitary only for real z; ', ...
          'all channels are retained; constant wide-band coupling and ', ...
          'no thresholds/absorption. Poles differ from eig(H).']);
    return
end
if ~(isnumeric(E)&&isscalar(E)&&isreal(E)&&isfinite(E))
    error('resonant_scattering:Energy', ...
        'The on-shell scattering energy must be real and finite.');
end
if strcmp(mode,'spectral')
    S=um_resonant_scattering_spectral(n,E,H,W,Gamma);
    return
end
S=um_resonant_scattering_eval(E,H,W,mode);
end
