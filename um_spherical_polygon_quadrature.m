function [K,areaNumerical]=um_spherical_polygon_quadrature(L,G,q)
%UM_SPHERICAL_POLYGON_QUADRATURE Positive radial-projection Gaussian rule.
%   K=sum_{triangles} int_T Y(x) Y(x)' dOmega, using real normalized
%   Condon--Shortley spherical harmonics in the same order as
%   SPHERICAL_SLEPIAN: m=0, then cos(m*phi), sin(m*phi) for each m>0;
%   within each (m,kind) the degrees ell run from m to L-1.
%
%   For a positively oriented spherical triangle with vertices a,b,c:
%     x(s,t) = (1-s)*a+s*(1-t)*b+s*t*c,       (s,t) in [0,1]^2
%     y(s,t) = x(s,t)/norm(x(s,t))
%     dOmega = s*det([a,b,c])/norm(x(s,t))^3 ds dt.
%   This is the radial projection of a planar simplex after a Duffy
%   square-to-simplex map. The Jacobian is ANALYTIC; the tensor
%   q-by-q Gauss-Legendre evaluation is an approximation.
%   It is NOT polynomial-exact in q and is NOT ARPIST's
%   rigorously bounded floating-point cubature implementation.
%   Each quadrature weight is positive so K is PSD up to roundoff.
%
%   [172] https://doi.org/10.1016/j.cam.2022.114822
%   [173] https://doi.org/10.1016/j.amc.2025.129335
validateattributes(q,{'numeric'}, ...
    {'scalar','real','integer','positive','finite','<=',96});
j=(1:q-1)';
off=j./sqrt(4*j.^2-1);
J=diag(off,1)+diag(off,-1);
[V,D]=eig(J);
[x,idx]=sort(real(diag(D)));
oneD=(x+1)/2;
w=real(V(1,idx).').^2; % normalized 1D weights on [0,1]
S=repmat(oneD,1,q);
T=repmat(oneD.',q,1);
s=S(:);t=T(:);
flatWeights=repmat(w,1,q).*repmat(w.',q,1).*S;
baseWeights=flatWeights(:);
n=L^2;
K=zeros(n);
areaNumerical=0;
for tri=1:size(G.triangles,1)
    v=G.vertices(G.triangles(tri,:),:);
    a=v(1,:)';b=v(2,:)';c=v(3,:)';
    det3=dot(a,cross(b,c));
    % 3-by-q² unnormalized chord points, with positive Jacobian.
    xyz=a*(1-s).'+b*(s.*(1-t)).'+c*(s.*t).';
    h=sqrt(sum(xyz.^2,1));
    xyz=xyz./h;
    weight=baseWeights.*(det3./h(:).^3);
    areaNumerical=areaNumerical+sum(weight);
    z=xyz(3,:);
    phi=atan2(xyz(2,:),xyz(1,:));
    B=zeros(q*q,n);
    Pl=cell(L,1);
    for ell=0:L-1,Pl{ell+1}=legendre(ell,z);end
    index=0;
    for m=0:L-1
        if m==0,kinds=0;else,kinds=[1 2];end
        for kind=kinds
            if kind==0,angular=ones(q*q,1)/sqrt(2*pi);
            elseif kind==1,angular=cos(m*phi(:))/sqrt(pi);
            else,angular=sin(m*phi(:))/sqrt(pi);
            end
            for ell=m:L-1
                f=sqrt((2*ell+1)/2*exp( ...
                    gammaln(ell-m+1)-gammaln(ell+m+1)));
                radial=f*reshape(Pl{ell+1}(m+1,:),[],1);
                index=index+1;
                B(:,index)=radial.*angular;
            end
        end
    end
    if index~=n
        error('um_spherical_polygon_quadrature:Internal', ...
            'Harmonic ordering dimension mismatch.');
    end
    BW=sqrt(weight).*B;
    K=K+BW'*BW;
end
end
