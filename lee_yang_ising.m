function [U,properties] = lee_yang_ising(n,couplings,pattern,mode)
%LEE_YANG_ISING Unitary Floquet-CMV matrix for Lee-Yang Ising zeros.
%   U=LEE_YANG_ISING(N) returns a sparse unitary matrix of order 2*N
%   whose eigenvalues are the Lee-Yang zeros for an N-site ferromagnetic
%   Ising chain with periodic spin boundary conditions. The default
%   dimensionless couplings p_j=J_j/(k_B*T) follow the Fibonacci word
%   a->ab, b->a, with p_a=2/3, p_b=1/100.
%
%   U=LEE_YANG_ISING(N,COUPLINGS,PATTERN) accepts:
%     'fibonacci'   COUPLINGS=[p_a,p_b] (the default)
%     'alternating' COUPLINGS=[p_a,p_b]
%     'uniform'     COUPLINGS=a positive scalar
%     'custom'      COUPLINGS=a positive N-vector [p_1,...,p_N].
%   A scalar is also accepted for either two-letter pattern.
%
%   All p_j must be positive. Set alpha(2*j-1)=exp(-2*p_j),
%   alpha(2*j)=0 (MATLAB one-based indexing). Form the periodic
%   2N-by-2N Floquet-CMV matrix F_{2N}(pi/2)=L*M, where L consists
%   of the even-numbered zero-based CMV cores and M of the odd cores,
%   including the pi/2-twisted cyclic boundary coupling. Unlike the
%   ordinary finite CMV truncation, every alpha_j lies in the disk.
%
%   U=LEE_YANG_ISING(...,'pencil') returns struct P with fields A=M,
%   B=L', U=L*M, satisfying B\A=U. 'spectral' returns the standard
%   descriptor with the exact secular characterization (no invented
%   closed-form density). Default output is the sparse matrix.
%
%   The unit-circle eigenvalues are precisely the zeros of
%     Z_N(z)=sum_{s_j=+/-1} exp(sum_j p_j*s_j*s_{j+1}) * z^(sum_j s_j),
%   with s_{N+1}=s_1. In particular z^N Z_N(z) is an even polynomial
%   of degree 2N. For Fibonacci couplings the thermodynamic zero
%   distribution can have Cantor-type gaps with gap labels.
%
%   Reference: D. Damanik, M. Embree and J. Fillman, Gap labels for
%   zeros of the partition function of the 1D Ising model via the
%   Schwartzman homomorphism, Indag. Math. 35 (2024), 813-836,
%   DOI 10.1016/j.indag.2023.05.004; Secs. 2-3, Example 5.2.
%   Background: A. Bunse-Gerstner and L. Elsner, Schur parameter pencils,
%   DOI 10.1016/0024-3795(91)90402-I.
%
%   Standard [U,PROPERTIES] signature preserves Anymatrix discovery.
properties={'square','unitary','normal','sparse','scalable','parameter-dependent'};
if nargin==0
    U=[];
    return
end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<3 || isempty(pattern),pattern='fibonacci';end
if ~(ischar(pattern)||isstring(pattern)) || ~isscalar(string(pattern))
    error('lee_yang_ising:Pattern','PATTERN must be a scalar string.');
end
pattern=lower(char(pattern));
if nargin<2 || isempty(couplings)
    if strcmp(pattern,'uniform'),couplings=2/3;
    else,couplings=[2/3 1/100];end
end
validateattributes(couplings,{'numeric'},{'vector','real','finite','positive'});
switch pattern
    case 'uniform'
        if numel(couplings)~=1
            error('lee_yang_ising:Couplings','Uniform pattern needs one coupling.');
        end
        p=repmat(double(couplings),n,1);
    case 'custom'
        if numel(couplings)~=n
            error('lee_yang_ising:Couplings','Custom pattern needs N couplings.');
        end
        p=double(couplings(:));
    case {'fibonacci','alternating'}
        if numel(couplings)==1,couplings=[couplings couplings];end
        if numel(couplings)~=2
            error('lee_yang_ising:Couplings','Two-letter pattern needs two couplings.');
        end
        if strcmp(pattern,'alternating')
            idx=mod((0:n-1)',2)+1;
        else
            % Finite prefix of S^k(a), S(a)=ab, S(b)=a.
            word='a';
            while numel(word)<n
                parts=cell(1,numel(word));
                for j=1:numel(word)
                    if word(j)=='a',parts{j}='ab';
                    else,parts{j}='a';end
                end
                word=[parts{:}]; %#ok<AGROW>
            end
            idx=1+double(word(1:n)=='b');
            idx=idx(:);
        end
        p=double(couplings(idx));
    otherwise
        error('lee_yang_ising:Pattern', ...
              'Supported patterns: fibonacci, alternating, uniform, custom.');
end
if nargin<4 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||isstring(mode)) || ~isscalar(string(mode))
    error('lee_yang_ising:Mode','MODE must be a scalar string.');
end
mode=lower(char(mode));
if strcmp(mode,'spectral')
    U=unitary_spectral_reference('finite_cmv',2*n);
    U.family='lee_yang_ising';
    U.reference='Damanik, Embree and Fillman, Indag. Math. 35 (2024), 813-836';
    U.notes=['Finite eigenvalues exactly coincide with the Lee-Yang zeros ' ...
        'of the stated Ising model. Fibonacci coupling yields gap-labelled ' ...
        'spectral features; no universal analytic trace density is asserted.'];
    return
end
if ~any(strcmp(mode,{'matrix','pencil'}))
    error('lee_yang_ising:Mode','Supported modes: matrix, pencil, spectral.');
end
N=2*n;
alpha=zeros(N,1);
alpha(1:2:end)=exp(-2*p);
rho=sqrt((1-alpha).*(1+alpha));
L=speye(N);
M=speye(N);
for k=1:2:N
    L(k:k+1,k:k+1)=[alpha(k),rho(k);rho(k),-alpha(k)];
end
for k=2:2:N-2
    M(k:k+1,k:k+1)=[alpha(k),rho(k);rho(k),-alpha(k)];
end
% The odd CMV core crossing the periodic boundary has twist exp(i*pi/2).
theta=pi/2;
M([1 N],[1 N])=[-alpha(N),exp(-1i*theta)*rho(N); ...
                   exp(1i*theta)*rho(N),alpha(N)];
U=L*M;
if strcmp(mode,'pencil')
    U=struct('A',M,'B',L','U',U);
end
end
