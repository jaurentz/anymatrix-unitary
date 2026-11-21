function [A,b,c,order]=um_polar_galerkin_coefficients(stages)
%UM_POLAR_GALERKIN_COEFFICIENTS Gauss--Legendre tableaux in Shen et al.
%   Stages=1,2,3 give one-, two-, and three-stage Gauss--Legendre
%   rules in Shen, Tran, Leok (2022), Table 1(A),(C),(D), with
%   classical collocation orders 2,4,6 respectively.
%
%   These are the tableaux USED WITHIN the variational polar
%   decomposition integrator, NOT a claim that a polar kinematic
%   stage alone is a symplectic/Lie--Poisson integrator.
%   Source: https://doi.org/10.3934/jcd.2022003
if ~(isnumeric(stages)&&isscalar(stages)&&isreal(stages) ...
        &&isfinite(stages)&&any(stages==[1 2 3]))
    error('um_polar_galerkin_coefficients:Stages', ...
        'Stages must be 1, 2 or 3.');
end
switch stages
    case 1
        A=.5;b=1;c=.5;order=2;
    case 2
        q=sqrt(3)/6;
        A=[.25 .25-q;.25+q .25];
        b=[.5;.5];
        c=[.5-q;.5+q];order=4;
    otherwise
        q=sqrt(15);
        A=[5/36 2/9-q/15 5/36-q/30; ...
           5/36+q/24 2/9 5/36-q/24; ...
           5/36+q/30 2/9+q/15 5/36];
        b=[5/18;4/9;5/18];
        c=[.5-q/10;.5;.5+q/10];order=6;
end
end
