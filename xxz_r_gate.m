function G=xxz_r_gate(u,gamma)
%XXZ_R_GATE Normalized unitary checked six-vertex/XXZ R-matrix.
%   Real rapidity u and anisotropy 0<gamma<pi, Delta=cos(gamma).
%   G is 4-by-4 in basis |00>,|01>,|10>,|11>.
%   R(u) has diagonal a=sinh(u+i*gamma), a, and central
%   [[sinh(u),i*sin(gamma)],[i*sin(gamma),sinh(u)]].
%   G=P*R(u)/sqrt(sinh(u)^2+sin(gamma)^2) with P the swap.
%   Therefore G'*G=I and G satisfies braid-form Yang-Baxter.
%   Internal helper, not an Anymatrix matrix generator.
%   Source: J.-S. Caux, The Bethe Ansatz, trigonometric R-matrix.
validateattributes(u,{'numeric'},{'scalar','real','finite'});
validateattributes(gamma,{'numeric'},{'scalar','real','finite','>',0,'<',pi});
s=sinh(u);
t=sin(gamma);
rho=hypot(s,t);
a=sinh(u+1i*gamma)/rho;
G=[a 0 0 0;0 1i*t/rho s/rho 0; ...
   0 s/rho 1i*t/rho 0;0 0 0 a];
end
