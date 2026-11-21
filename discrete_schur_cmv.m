function [U,properties]=discrete_schur_cmv(n,k,h,nodes,weights,output)
%DISCRETE_SCHUR_CMV Finite positive-measure steps for a CMV Schur-flow benchmark.
%   U=DISCRETE_SCHUR_CMV(N,K,H,NODES,WEIGHTS,OUTPUT) applies K integer
%   positive Christoffel/Laurent measure steps to a fixed N-node circle
%   spectrum. For 0<H<1/2 and |z|=1, g_h(z)=1+H*(z+conj(z))>0 and
%
%        w_j(K) = w_j(0)*g_h(z_j)^K / sum_l w_l(0)*g_h(z_l)^K.
%
%   The support never changes. INVERSE_OPUC constructs an exactly
%   unitary finite CMV matrix with these cyclic e1 weights. Its finite
%   TRACE spectral measure instead assigns equal mass 1/N to each node.
%   All powers and the terminal phase are fixed by the support.
%
%   K=0 returns the original finite inverse-spectral matrix. The
%   continuous-time limit K*H -> T, H -> 0 yields weights proportional
%   to w_j(0)*exp(2*T*real(z_j)). The discrete weight deformation is a
%   mathematically exact positive finite-measure benchmark; it is NOT
%   claimed to be the Mukaihira--Nakamura integrable time discretization.
%
%   OUTPUT='matrix' (default), 'operator', 'parameters' or 'spectral'.
%   Operator mode uses the O(N)-storage CMV core product. Large K can
%   cause positive weights to underflow in floating point; in that case
%   inverse spectral reconstruction is rejected rather than altering
%   the support or claiming high accuracy.
%
%   Mukaihira--Nakamura, JCAM 139 (2002), 75--94:
%   https://doi.org/10.1016/S0377-0427(01)00388-0
%   Cantero et al., Adv. Math. (2016):
%   https://doi.org/10.1016/j.aim.2016.03.042
%   Simon, OPUC (2005): https://doi.org/10.1090/coll054.1
properties={'square','unitary','normal','scalable', ...
    'parameter-dependent','eigenvalues'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(k),k=0;end
if nargin<3||isempty(h),h=.1;end
validateattributes(k,{'numeric'},{'scalar','real','integer','nonnegative','finite'});
validateattributes(h,{'numeric'},{'scalar','real','finite','>',0,'<',.5});
if nargin<4||isempty(nodes)
    nodes=exp(2i*pi*((0:n-1)'+.17)/n);
end
if nargin<5||isempty(weights),weights=ones(n,1)/n;end
if nargin<6||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('discrete_schur_cmv:Mode','OUTPUT must be a scalar string.');
end
output=lower(char(output));
if ~any(strcmp(output,{'matrix','operator','parameters','spectral'}))
    error('discrete_schur_cmv:Mode','Invalid output mode.');
end
old=inverse_opuc(n,nodes,weights,'cmv','spectral');
z=old.finite.eigenvalues();
w=weights(:);
w=w/sum(w);
logf=log1p(2*h*real(z));
logw=log(w)+k*logf;
logw=logw-max(logw);
wk=exp(logw);
wk=wk/sum(wk);
if any(wk<=0)||any(~isfinite(wk))
    error('discrete_schur_cmv:WeightUnderflow', ...
        'Strictly positive finite weights underflowed; decrease K.');
end
switch output
    case 'spectral'
        U=inverse_opuc(n,z,wk,'cmv','spectral');
        U.family='discrete_schur_cmv';
        U.discrete_schur=struct('step',k,'h',h,'initial_weights',w, ...
            'weights',wk,'factors',1+2*h*real(z), ...
            'time_limit',k*h,'kind','positive_finite_measure_steps');
        U.reference='Mukaihira--Nakamura (2002); Cantero et al. (2016)';
        U.notes=[U.notes ' Each finite CMV spectrum is purely atomic. ' ...
            'This positive-measure step is not identified with the ' ...
            'published integrable Mukaihira--Nakamura discretization.'];
    case 'parameters'
        U=inverse_opuc(n,z,wk,'cmv','parameters');
        U.step=k;
        U.h=h;
        U.initial_weights=w;
        U.factors=1+2*h*real(z);
        U.time_limit=k*h;
    otherwise
        U=inverse_opuc(n,z,wk,'cmv',output);
end
end
