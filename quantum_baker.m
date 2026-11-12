function [U,properties] = quantum_baker(n,variant)
%QUANTUM_BAKER Finite quantized baker map (Balazs-Voros or Saraceno).
%   U=QUANTUM_BAKER(N) is F_N'*diag(F_(N/2),F_(N/2)), N positive even,
%   where F_N(j,k)=exp(-2*pi*1i*j*k/N)/sqrt(N).
%   U=QUANTUM_BAKER(N,'saraceno') uses half-integer Fourier phases in
%   both full and half-size blocks. Operator and boundary conventions
%   must be recorded when comparing eigenphase statistics.
%   References: Balazs-Voros, Ann. Phys. 190 (1989), 1-31;
%   Saraceno, Ann. Phys. 199 (1990), 37-60.
properties={'square','unitary','normal','complex','scalable','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(n,{'numeric'},{'scalar','integer','positive','finite'});
if mod(n,2),error('quantum_baker:EvenOrder','N must be even.');end
if nargin<2 || isempty(variant),variant='balazs_voros';end
variant=lower(char(variant));
if ~any(strcmp(variant,{'balazs_voros','saraceno'}))
    error('quantum_baker:Variant','Unknown quantization convention.');
end
shift=double(strcmp(variant,'saraceno'))/2;
F=local_fourier(n,shift);
G=local_fourier(n/2,shift);
U=F'*blkdiag(G,G);
end

function F=local_fourier(n,shift)
q=(0:n-1)+shift;
F=exp(-2*pi*1i*(q(:)*q(:).')/n)/sqrt(n);
end
