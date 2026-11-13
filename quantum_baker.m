function [U,properties] = quantum_baker(n,variant,output)
%QUANTUM_BAKER Finite quantized baker map (Balazs-Voros or Saraceno).
%   U=QUANTUM_BAKER(N) is F_N'*diag(F_(N/2),F_(N/2)), N positive even,
%   where F_N(j,k)=exp(-2*pi*1i*j*k/N)/sqrt(N).
%   U=QUANTUM_BAKER(N,'saraceno') uses half-integer Fourier phases in
%   both full and half-size blocks. Operator and boundary conventions
%   must be recorded when comparing eigenphase statistics.
%   References: Balazs-Voros, Ann. Phys. 190 (1989), 1-31;
%   Saraceno, Ann. Phys. 199 (1990), 37-60.
%   Verified reference identifiers (see docs/MATRIX_ENCYCLOPEDIA.md):
%   [20] https://doi.org/10.1016/0003-4916(89)90259-5
%   [21] https://doi.org/10.1016/0003-4916(90)90367-W
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
% Optional analytic spectral-reference mode.
if nargin>=3 && (ischar(output)||isstring(output)) && strcmpi(char(output),'spectral')
    U=unitary_spectral_reference('quantum_baker',n,variant);
    return
end
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if mod(n,2),error('quantum_baker:EvenOrder','N must be even.');end
if nargin<2 || isempty(variant),variant='balazs_voros';end
variant=lower(char(variant));
if ~any(strcmp(variant,{'balazs_voros','saraceno'}))
    error('quantum_baker:Variant','Unknown quantization convention.');
end
shift=double(strcmp(variant,'saraceno'))/2;
if nargin<3,output='matrix';end
op=@(X,flag)um_fft_apply(X,flag,'baker',n,shift);
if um_is_operator(output),U=op;else,U=op(eye(n),'notransp');end
end
