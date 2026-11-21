function tests=test_hankel_nonunitarity
%TEST_HANKEL_NONUNITARITY Guard against false finite DHT unitarity.
%   The usual Bessel-zero quasi-DHT is generally NEARLY orthogonal, not
%   exactly orthogonal at finite matrix order. This is deliberately a
%   negative test, not a generator registered as unitary in Anymatrix.
%
%   Baddour--Chouinard, JOSA A 32 (2015), 611--622, Eq. (44):
%     T(j,k)=2*J_nu(z_j*z_k/z_{N+1})/
%       (z_{N+1}*J_{nu+1}(z_j)*J_{nu+1}(z_k)).
%   Their underlying finite Bessel orthogonality identity is explicitly
%   identified as asymptotic (see Baddour 2019, sec. 4.4).
%   https://doi.org/10.1364/JOSAA.32.000611
%   https://www.intechopen.com/chapters/65719
tests=functiontests(localfunctions);
end

function test_bessel_zero_truncations_are_not_exact_involutions(tc)
% All three positive integer orders are deliberately distinct from the
% elementary half-integer Bessel/DST reduction.
for nu=[0 1 2]
    for n=[3 5 8]
        T=raw_bessel_zero_dht(n,nu);
        verifyEqual(tc,size(T),[n n]);
        verifyLessThan(tc,norm(T-T.','fro'),1e-13*n);
        gram=norm(T.'*T-eye(n),'fro');
        invol=norm(T*T-eye(n),'fro');
        % Not numerical roundoff: true finite truncation discrepancies.
        verifyGreaterThan(tc,gram,1e-8);
        verifyGreaterThan(tc,invol,1e-8);
        % These are close-to-involutory, but NOT orthogonal matrices.
        verifyLessThan(tc,gram,1e-2);
        verifyLessThan(tc,invol,1e-2);
    end
end
end

function test_half_integer_exact_dst_i_specialization(tc)
% J_{1/2}(x)=sqrt(2/(pi*x))*sin(x), so its first positive zeros
% are k*pi and the same Bessel-zero normalization is exactly DST-I.
% This is an existing orthogonal transform, not a new Bessel family.
nu=.5;
for n=[3 5 8]
    a=pi*(1:n)';
    c=(n+1)*pi;
    d=abs(besselj(nu+1,a));
    T=(2/c)*besselj(nu,(a*a.')/c)./(d*d.');
    j=(1:n)';
    DST=sqrt(2/(n+1))*sin(pi*(j*j.')/(n+1));
    verifyLessThan(tc,norm(T-DST,'fro'),2e-12);
    verifyLessThan(tc,norm(T'*T-eye(n),'fro'),2e-12);
end
end

function test_small_order_reproducible_defect(tc)
T=raw_bessel_zero_dht(5,0);
% Independent scipy.special.jn_zeros / jv evaluation gives
% ||T^T T-I||_F = 3.6356e-6 for this explicitly stated normalization.
verifyLessThan(tc,abs(norm(T.'*T-eye(5),'fro') ...
    -3.6356e-6),1e-9);
end

function T=raw_bessel_zero_dht(n,nu)
z=positive_bessel_zeros(nu,n+1);
a=z(1:n);
b=z(n+1);
d=abs(besselj(nu+1,a));
T=(2/b)*besselj(nu,(a*a.')/b)./(d*d.');
end

function roots=positive_bessel_zeros(nu,count)
% Independent elementary sign-bracketing of integer-order Bessel roots.
% No Bessel zero or special-function toolbox is needed.
roots=zeros(count,1);
h=pi/8;
left=1e-8;
fleft=besselj(nu,left);
found=0;
for k=1:5000
    right=left+h;
    fright=besselj(nu,right);
    if fleft*fright<0
        found=found+1;
        roots(found)=fzero(@(x)besselj(nu,x),[left,right]);
        if found==count,return;end
    end
    left=right;
    fleft=fright;
end
error('test_hankel_nonunitarity:Roots', ...
    'Could not bracket the requested positive Bessel zeros.');
end
