function [U,properties]=fibonacci_braid(n,charge,word,mode,chirality)
%FIBONACCI_BRAID Fixed-charge Fibonacci-anyon braid matrices.
%   U=FIBONACCI_BRAID(N,CHARGE,WORD,MODE,CHIRALITY) acts on the
%   fusion-path space of N Fibonacci tau anyons, N=3,...,20.
%   CHARGE='tau' (default) or 'vacuum' is the conserved total charge.
%   Its dimension is Fibonacci(N) for tau and Fibonacci(N-1)
%   for vacuum (F1=F2=1). Basis vectors are left-associated
%   fusion paths [a1,...,aN], with a1=tau, aN=CHARGE,
%   a_j in {vacuum,tau}, and fusion rules tau*tau=vacuum+tau,
%   vacuum*tau=tau. Paths are in lexicographic order, vacuum=0.
%
%   WORD is a signed integer vector of braid indices: +i exchanges
%   anyons i and i+1 counterclockwise, -i is its inverse, and the
%   first word entry acts first. [] is the identity braid. Default
%   WORD=[1 2 -1]. MODE is 'matrix' (default, dimension<=512),
%   'operator' (sparse factor action/adjoint), 'spectral' (analytic
%   finite trace spectrum for identity and one generator), or
%   'generators' (fusion paths and sparse sigma_i matrices).
%   CHIRALITY=+1 (default) uses
%       R(vacuum)=exp(-4*pi*i/5), R(tau)=exp(3*pi*i/5).
%   CHIRALITY=-1 complex conjugates these R phases. The nontrivial
%   F-symbol is a real symmetric involution with entries
%       [1/phi, 1/sqrt(phi);1/sqrt(phi),-1/phi],
%   phi=(1+sqrt(5))/2. This is a specified fusion-tree gauge;
%   changing basis gauge conjugates all braid matrices, but changing
%   chirality conjugates physical phases and changes eigenphases.
%
%   This is an actual non-Abelian braid representation, NOT a
%   tensor-product qubit gate circuit or a generic unitary sample.
%   S_i obey S_i*S_(i+1)*S_i=S_(i+1)*S_i*S_(i+1) and commute
%   for distant indices. Each generator is sparse with O(dim)
%   nonzeros, and the operator action applies word factors in order.
%
%   S. Trebst et al., Prog. Theor. Phys. Suppl. 176 (2008)
%   384-407 [130]. https://doi.org/10.1143/PTPS.176.384
%   C. Carnahan, D. Zeuch, N. E. Bonesteel, Phys. Rev. A
%   93 (2016) 052328 [131].
%   https://doi.org/10.1103/PhysRevA.93.052328
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','real','integer','finite','>=',3,'<=',20});
if nargin<2 || isempty(charge),charge='tau';end
if nargin<3,word=[1 2 -1];end
if nargin<4 || isempty(mode),mode='matrix';end
if nargin<5 || isempty(chirality),chirality=1;end
if ~(ischar(charge)||(isstring(charge)&&isscalar(charge)))
    error('fibonacci_braid:Charge','Charge must be text.');
end
charge=lower(strtrim(char(charge)));
switch charge
    case {'tau','t'}
        total=1;charge='tau';
    case {'vacuum','1','one'}
        total=0;charge='vacuum';
    otherwise
        error('fibonacci_braid:Charge','Charge must be tau or vacuum.');
end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('fibonacci_braid:Mode','Mode must be text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','operator','spectral','generators'}))
    error('fibonacci_braid:Mode','Unknown output mode.');
end
if isempty(word)
    word=zeros(1,0);
elseif ~(isnumeric(word)&&isreal(word)&&isvector(word) ...
         && all(isfinite(word)) && all(word==round(word)) ...
         && all(word~=0) && all(abs(word)<=n-1))
    error('fibonacci_braid:Word','Word must be signed braid indices 1,...,N-1.');
end
word=double(word(:).');
if ~(isnumeric(chirality)&&isscalar(chirality)&&isreal(chirality) ...
       && isfinite(chirality) && any(chirality==[-1 1]))
    error('fibonacci_braid:Chirality','Chirality must be +1 or -1.');
end
paths=fusion_paths(n,total);
d=size(paths,1);
phi=(1+sqrt(5))/2;
F=[1/phi,1/sqrt(phi);1/sqrt(phi),-1/phi];
r=[exp(-chirality*4i*pi/5);exp(chirality*3i*pi/5)];
if strcmp(mode,'spectral')
    U=um_fibonacci_braid_spectral(n,charge,word,chirality,paths,r);
    return
end
pow=2.^(n-1:-1:0)';
codes=1+paths*pow;
lookup=sparse(codes,ones(d,1),1:d,2^n,1);
FRF=F*diag(r)*F;
sig=cell(n-1,1);
for k=1:n-1
    sig{k}=local_sigma(paths,lookup,pow,k,r,FRF);
end
if strcmp(mode,'generators')
    U=struct('braid_generators',{sig},'fusion_paths',paths, ...
        'dimension',d,'anyon_count',n,'charge',charge, ...
        'chirality',chirality,'F',F,'R',r);
    return
end
O=@(X,flag)um_fibonacci_braid_apply(X,flag,d,sig,word);
if strcmp(mode,'operator')
    U=O;
else
    if d>512
        error('fibonacci_braid:ExplicitOrder', ...
            'Dense output limited to fusion-space dimension <=512.');
    end
    U=O(eye(d),'notransp');
end
end

function paths=fusion_paths(n,total)
paths=1; % a1=tau; vacuum 0, tau 1
for step=2:n
    next=zeros(0,step);
    for j=1:size(paths,1)
        p=paths(j,:);
        if p(end)==0
            next(end+1,:)=[p 1]; %#ok<AGROW>
        else
            next(end+1,:)=[p 0]; %#ok<AGROW>
            next(end+1,:)=[p 1]; %#ok<AGROW>
        end
    end
    paths=next;
end
paths=paths(paths(:,end)==total,:);
end

function G=local_sigma(paths,lookup,pow,k,r,FRF)
d=size(paths,1);
rows=zeros(2*d,1);cols=zeros(2*d,1);
values=complex(zeros(2*d,1));count=0;
for col=1:d
    p=paths(col,:);
    if k==1
        count=count+1;
        rows(count)=col;cols(count)=col;values(count)=r(p(2)+1);
    else
        left=p(k-1);mid=p(k);right=p(k+1);
        if left==0
            count=count+1;
            rows(count)=col;cols(count)=col;values(count)=r(right+1);
        elseif right==0
            count=count+1;
            rows(count)=col;cols(count)=col;values(count)=r(2);
        else
            for outcome=0:1
                p(k)=outcome;
                row=full(lookup(1+p*pow));
                if row<1
                    error('fibonacci_braid:Fusion','Invalid recoupling path.');
                end
                count=count+1;
                rows(count)=row;cols(count)=col;
                values(count)=FRF(outcome+1,mid+1);
            end
        end
    end
end
G=sparse(rows(1:count),cols(1:count),values(1:count),d,d);
end
