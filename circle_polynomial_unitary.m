function [U,properties]=circle_polynomial_unitary(n,coeffs,weights,representation,output)
%CIRCLE_POLYNOMIAL_UNITARY Unitary companion-like circle polynomial model.
%   U=CIRCLE_POLYNOMIAL_UNITARY(N,COEFFS,WEIGHTS,REP,OUTPUT)
%   constructs a strict N-by-N unitary Hessenberg or finite CMV
%   realization of a degree-N polynomial having N DISTINCT roots
%   on the unit circle. COEFFS are descending monomial coefficients;
%   the polynomial is normalized to monic form before realization.
%   REP='hessenberg' (default) or 'cmv'; WEIGHTS are positive cyclic
%   spectral weights (default equal) and need not be trace weights.
%   This calls INVERSE_OPUC for the actual Verblunsky/core products.
%   Neither the ordinary monomial companion nor an arbitrary
%   self-inversive polynomial is automatically unitary.
%   The input must be approximately self-inversive AND all its
%   numerically computed, simple roots must lie on the unit circle.
%   Ill-conditioned or near-multiple roots can be rejected.
%   'operator' uses existing OPUC core actions and adjoint.
%   'parameters' includes the normalized monic polynomial, roots,
%   reflection coefficients and closure phase.
%   'pencil' uses the existing strict finite block-CMV pencil, with
%   the same eigenvalues but not necessarily the Hessenberg matrix.
%   'spectral' reports the constructed finite trace and e1 measures.
%   ROOTS is a numerical polynomial solver; no exact phase formula
%   is claimed from symbolic polynomial data.
%   [5] https://doi.org/10.1090/coll054.1
%   [5] https://doi.org/10.1090/coll/054.2
%   [9] https://doi.org/10.1016/j.cam.2014.11.002
properties={'square','unitary','normal','hessenberg','scalable', ...
    'parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if n>64
    error('circle_polynomial_unitary:Order', ...
        'Polynomial-root-to-CMV inversion capped at n=64 for reliability.');
end
if nargin<2||isempty(coeffs)
    coeffs=[1,zeros(1,n-1),-exp(1i*n*.11)];
end
if ~isnumeric(coeffs)||~isvector(coeffs)||numel(coeffs)~=n+1 ...
        ||any(~isfinite(coeffs(:)))||coeffs(1)==0
    error('circle_polynomial_unitary:Polynomial', ...
        'COEFFS must be N+1 finite descending coefficients of degree N.');
end
coeffs=coeffs(:).'/coeffs(1);
if abs(abs(coeffs(end))-1)>2e-9 || ...
        norm(coeffs-coeffs(end)*conj(fliplr(coeffs)),inf)> ...
        3e-9*max(1,norm(coeffs,inf))
    error('circle_polynomial_unitary:SelfInversive', ...
        'A circle-root polynomial must be self-inversive.');
end
z=roots(coeffs);
if any(~isfinite(z))||any(abs(abs(z)-1)>5e-8) || ...
        max(abs(polyval(coeffs,z)))> ...
        5e-8*max(1,norm(coeffs,1))
    error('circle_polynomial_unitary:CircleRoots', ...
        'All polynomial roots must lie on the unit circle.');
end
z=z./abs(z);  % remove small root solver radial errors
for j=2:n
    if min(abs(z(j)-z(1:j-1)))<=128*eps*max(1,n)
        error('circle_polynomial_unitary:Multiplicity', ...
            'Distinct simple circle roots are required.');
    end
end
if nargin<3||isempty(weights),weights=ones(n,1)/n;end
if nargin<4||isempty(representation),representation='hessenberg';end
if ~(ischar(representation)||(isstring(representation)&&isscalar(representation)))
    error('circle_polynomial_unitary:Representation', ...
        'REP must be hessenberg or cmv.');
end
representation=lower(char(representation));
if ~any(strcmp(representation,{'hessenberg','cmv'}))
    error('circle_polynomial_unitary:Representation', ...
        'REP must be hessenberg or cmv.');
end
if nargin<5||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('circle_polynomial_unitary:Mode','Invalid output mode.');
end
switch lower(char(output))
    case 'pencil'
        info=inverse_opuc(n,z,weights,'cmv','parameters');
        U=block_cmv_pencil(n,1,reshape(info.alpha,1,1,[]), ...
            info.beta,'pencil');
        U.polynomial=coeffs;
        U.roots=z;
    case 'parameters'
        U=inverse_opuc(n,z,weights,representation,'parameters');
        U.polynomial=coeffs;
    case 'spectral'
        U=inverse_opuc(n,z,weights,representation,'spectral');
        U.polynomial=coeffs;
        U.notes=[U.notes ' Phases derived by numerical roots of the ' ...
            'supplied polynomial; no exact symbolic root formula claimed.'];
    case {'matrix','explicit','operator','implicit'}
        U=inverse_opuc(n,z,weights,representation,output);
    otherwise
        error('circle_polynomial_unitary:Mode','Unknown output mode.');
end
end
