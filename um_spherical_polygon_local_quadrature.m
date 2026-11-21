function [K,areaNumerical,diag]=um_spherical_polygon_local_quadrature(L,G,tol)
%UM_SPHERICAL_POLYGON_LOCAL_QUADRATURE Locally adapted Gaussian polygon rule.
%   [K,AREA,D]=UM_SPHERICAL_POLYGON_LOCAL_QUADRATURE(L,G,TOL)
%   integrates a bandlimited concentration matrix over each of the
%   spherical triangles returned by UM_SPHERICAL_POLYGON_GEOMETRY.
%   The source integration kernel is UM_SPHERICAL_POLYGON_QUADRATURE.
%   Each triangle independently doubles its tensor-Gauss order
%   q=max(8,L+3), 2*q, ... up to 64. Accepted local tolerances
%   are apportioned so the SUM of triangle indicators obeys
%   the requested global absolute Frobenius scale. This is useful
%   when only a small number of polygon triangles need refinement.
%
%   The numerical error INDICATOR for triangle i is the Frobenius
%   difference ||Ki(q)-Ki(q/2)||F; the independent geometry check
%   compares cubature area against the EXACT spherical excess
%      Ai = 2 atan2(|det(a,b,c)|, 1+a.b+b.c+c.a).
%   By the matrix triangle inequality, the sum of measured local
%   differences bounds the difference between the assembled
%   high-order and assembled previous-order rules. IT DOES NOT
%   bound error to the mathematical integral. No certified claim.
%
%   Return diagnostic struct even when no local order reaches
%   tolerance (K,AREA are then last approximations); caller
%   decides whether to throw or expose achieved estimates.
%   Symmetry and spectral PSD/contraction must be checked by caller.
%   Computing each K uses positive-weight radial-projection
%   Gauss-Legendre quadrature, preserving the mathematical
%   Condon--Shortley order from spherical_polygon_slepian.
%
%   Background: Li--Jiao, JCAM 2023 (radial projection);
%   Sommariva--Vianello, AMC 2025 (spherical polygon cubature).
nTri=size(G.triangles,1);
N=L^2;
K=zeros(N);
areaNumerical=0;
qStart=max(8,L+3);
qMax=64;
targetTrace=N*G.area/(4*pi);
globalScale=max(1,targetTrace);
% Conservative per-triangle allocation of the GLOBAL requested
% norm tolerance; the area comparison uses the same scale.
localGoal=tol*globalScale/max(1,nTri);
localOrders=zeros(nTri,1);
localLevels=zeros(nTri,1);
localDifferences=inf(nTri,1);
localAreaErrors=inf(nTri,1);
localExactAreas=zeros(nTri,1);
localConverged=false(nTri,1);
for tri=1:nTri
    idx=G.triangles(tri,:);
    a=G.vertices(idx(1),:);
    b=G.vertices(idx(2),:);
    c=G.vertices(idx(3),:);
    det3=dot(a,cross(b,c));
    den=1+dot(a,b)+dot(b,c)+dot(c,a);
    localExactAreas(tri)=2*atan2(abs(det3),den);
    sub=G;
    sub.triangles=idx;
    q=qStart;
    [previous,~]=um_spherical_polygon_quadrature(L,sub,q);
    localLevels(tri)=1;
    if 2*q>qMax
        error('um_spherical_polygon_local_quadrature:Internal', ...
            'Initial order leaves no room for refinement.');
    end
    while 2*q<=qMax
        q=2*q;
        [current,aNum]=um_spherical_polygon_quadrature(L,sub,q);
        localLevels(tri)=localLevels(tri)+1;
        localDifferences(tri)=norm(current-previous,'fro');
        localAreaErrors(tri)=abs(aNum-localExactAreas(tri));
        previous=current;
        if max(localDifferences(tri),localAreaErrors(tri))<=localGoal
            localConverged(tri)=true;
            break
        end
    end
    K=K+previous;
    areaNumerical=areaNumerical+aNum;
    localOrders(tri)=q;
end
matrixIndicator=sum(localDifferences);
areaError=abs(areaNumerical-G.area);
areaTriangleSumError=abs(sum(localExactAreas)-G.area);
traceError=abs(trace(K)-targetTrace);
converged=all(localConverged) && ...
    max([matrixIndicator,areaError,traceError])<=tol*globalScale;
diag=struct('triangle_orders',localOrders, ...
    'triangle_levels',localLevels, ...
    'triangle_frobenius_differences',localDifferences, ...
    'triangle_area_errors',localAreaErrors, ...
    'triangle_exact_areas',localExactAreas, ...
    'triangle_converged',localConverged, ...
    'total_matrix_difference_upper_sum',matrixIndicator, ...
    'total_area_error',areaError, ...
    'sum_exact_triangle_areas_error',areaTriangleSumError, ...
    'global_addition_trace_error',traceError, ...
    'global_scale',globalScale, ...
    'local_tolerance',localGoal, ...
    'requested_tolerance',tol, ...
    'converged',converged, ...
    'certified',false, ...
    'interpretation',['The sum of local rule differences is an ', ...
      'upper bound only on the difference between the TWO ', ...
      'assembled numerical Gaussian rules; it is NOT a bound ', ...
      'on error against the mathematical concentration integral.']);
end
