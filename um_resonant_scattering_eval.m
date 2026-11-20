function X=um_resonant_scattering_eval(z,H,W,mode)
%UM_RESONANT_SCATTERING_EVAL On-shell or meromorphic wide-band scattering.
%   S=I-2iW'*A(z)^(-1)W with A=z*I-H+iWW'.
%   For real z, Q=-iS'*dS/dz=2(A\W)'(A\W) is Hermitian PSD.
%   'matrix' and 'derivative' accept complex z away from poles;
%   'time_delay' requires real z and is not a holomorphic continuation.
if ~(isnumeric(z)&&isscalar(z)&&all(isfinite(z)))
    error('um_resonant_scattering_eval:Energy','Energy must be finite scalar.');
end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('um_resonant_scattering_eval:Mode','Mode must be text.');
end
mode=lower(char(mode));
if strcmp(mode,'time_delay') && ~isreal(z)
    error('um_resonant_scattering_eval:Energy', ...
        'The Hermitian Wigner-Smith delay matrix requires real energy.');
end
m=size(H,1);
A=z*eye(m)-H+1i*(W*W');
Y=A\W;
switch mode
    case 'matrix'
        X=eye(size(W,2))-2i*(W'*Y);
    case 'derivative'
        X=2i*W'*(A\Y);
    case 'time_delay'
        X=2*(Y'*Y);
        X=(X+X')/2; % Remove only floating-point asymmetry.
    otherwise
        error('um_resonant_scattering_eval:Mode','Unknown evaluation mode.');
end
end
