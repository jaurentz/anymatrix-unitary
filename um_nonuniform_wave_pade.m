function coeff=um_nonuniform_wave_pade(method)
%UM_NONUNIFORM_WAVE_PADE Diagonal Gauss-Legendre Pade numerator coefficients.
%   R_s(z)=P_s(z)/P_s(-z), with s stages and order 2s.
switch lower(char(method))
    case 'midpoint',coeff=[1 .5];
    case 'gauss2',coeff=[1 .5 1/12];
    case 'gauss3',coeff=[1 .5 1/10 1/120];
    otherwise
        error('um_nonuniform_wave_pade:Method','Unknown Gauss method.');
end
end
