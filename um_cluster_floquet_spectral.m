function S=um_cluster_floquet_spectral(n,theta,boundary,kick)
%UM_CLUSTER_FLOQUET_SPECTRAL Finite trace spectra at ideal pulse.
%   Open cluster: each set of n-2 independent K_j eigenvalues carries
%   (+exp(-i sum theta_j k_j), -exp(-i sum theta_j k_j)), each twice.
%   Periodic even cluster: C=prod CZ, CK_jC=X_j, and
%   P=prod_j K_j, so lambda=prod_j k_j*exp(-i sum theta_j k_j).
%   These are EXACT finite sample spectra, not ensemble laws.
N=2^n;
blank=struct('kind','unknown','density',[],'support',[], ...
    'atoms',[],'moments',[]);
S=struct('family','cluster_floquet','dimension',N, ...
    'normalization','probability (density with respect to dtheta)', ...
    'finite',struct('kind','unknown','eigenvalues',[], ...
                    'moments',[],'atoms',[]), ...
    'limit',blank,'ensemble',blank,'underlying',blank, ...
    'reference','Potter, Morimoto and Vishwanath (2016); Else and Nayak (2016)', ...
    'notes',['Finite TRACE spectrum of a Z2 x Z2 symmetric cluster ', ...
       'Floquet drive; exact only for ideal KICK=0. ', ...
       'Boundary-local pi-edge pairing is not implied for perturbed ', ...
       'or periodic drives; no CUE spectral law is claimed.']);
if kick~=0,return;end
if strcmp(boundary,'open')
    active=theta(2:n-1);
    S.finite.eigenvalues=@()open_eigenvalues(active);
    S.finite.moments=@(k)open_moments(active,k);
    S.notes=['Exact finite multiset: for each of 2^(n-2) cluster ', ...
        'stabilizer sign choices, phases +/-exp(-i sum theta_j k_j), ', ...
        'both twice repeated. Strict boundary pi pairing and zero ', ...
        'odd normalized trace moments hold for all bulk disorder.'];
else
    S.finite.eigenvalues=@()periodic_eigenvalues(theta);
    S.finite.moments=@(k)periodic_moments(theta,k);
    S.notes=['Even periodic cluster: exact finite multiset ', ...
        'lambda=(prod_j k_j)*exp(-i sum theta_j k_j), ', ...
        'one per stabilizer sign configuration. No boundary modes.'];
end
S.finite.kind='exact';
end

function z=open_eigenvalues(theta)
base=sign_eigenvalues(theta,false);
z=[base;base;-base;-base];
end

function z=periodic_eigenvalues(theta)
z=sign_eigenvalues(theta,true);
end

function z=sign_eigenvalues(theta,withProduct)
m=numel(theta);N=2^m;
bits=uint64((0:N-1)');
phase=zeros(N,1);prodSign=ones(N,1);
for j=1:m
    signj=1-2*double(bitget(bits,j));
    phase=phase+theta(j)*signj;
    if withProduct,prodSign=prodSign.*signj;end
end
z=prodSign.*exp(-1i*phase);
end

function moment=open_moments(theta,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
sz=size(k);kk=reshape(k,1,[]);
phase=prod(cos(bsxfun(@times,theta(:),kk)),1);
phase(mod(abs(kk),2)==1)=0;
moment=reshape(phase,sz);
end

function moment=periodic_moments(theta,k)
validateattributes(k,{'numeric'},{'real','integer','finite'});
sz=size(k);kk=reshape(k,1,[]);
c=prod(cos(bsxfun(@times,theta(:),kk)),1);
s=prod(-1i*sin(bsxfun(@times,theta(:),kk)),1);
c(mod(abs(kk),2)==1)=s(mod(abs(kk),2)==1);
moment=reshape(c,sz);
end
