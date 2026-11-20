function [U,properties]=paraunitary_fir(L,d,directions,U0,omega,output)
%PARAUNITARY_FIR Finite periodic M-channel paraunitary FIR filter bank.
%   U=PARAUNITARY_FIR(L,D,DIRECTIONS,U0,OMEGA,OUTPUT) builds a D*L
%   periodic block-circulant unitary analysis operator with L temporal
%   polyphase blocks of D channels. Column j of DIRECTIONS specifies
%   a normalized direction v_j of a degree-one FIR lattice section
%       S_j(z)=I+(z^(-1)-1)*v_j*v_j',  |z|=1.
%   Define E(z)=U0*S_M(z)*...*S_1(z), with U0 unitary.
%   The finite periodic time action applies a rank-one channel
%   projection at time t with a one-block circular delay x(t-1).
%   This is a genuine unitary D*L finite transform, not merely a
%   D-by-D unitary frequency response. Default D=2, L=8 and two
%   noncommuting real section directions.
%
%   OUTPUT='matrix': explicit D*L operator (guarded for large orders).
%   OUTPUT='operator': O(D*L*M) fast periodic action and adjoint.
%   OUTPUT='frequency': small E(exp(1i*OMEGA)), not a finite operator.
%   OUTPUT='polyphase': D-by-D-by-(M+1) Laurent coefficient array;
%       taps(:,:,r+1) is the coefficient of z^(-r), winding=-M.
%   OUTPUT='spectral': exact finite trace measure from the L Fourier
%       frequency blocks E(exp(2*pi*1i*k/L)), computed via D-by-D
%       eigensolves without constructing the full finite matrix.
%   The winding number of det(E(exp(1i*OMEGA))) is -M, and each
%   section is paraunitary on the entire unit circle.
%   This is distinct from BLASCHKE_POTAPOV (one rational response)
%   and DAUBECHIES_WAVELET (particular critically sampled QMF).
%   [15] https://doi.org/10.1137/1.9781611970104
%   [120] https://www.informit.com/store/multirate-systems-and-filter-banks-9780136057185
properties={'square','unitary','normal','complex','scalable', ...
    'eigenvalues','parameter-dependent'};
if nargin==0,U=[];return;end
validateattributes(L,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<2||isempty(d),d=2;end
validateattributes(d,{'numeric'},{'scalar','integer','positive','finite'});
if nargin<3||isequal(size(directions),[0 0])
    directions=zeros(d,2);
    directions(1,1)=1;
    directions(1,2)=1;
    directions(min(2,d),2)=directions(min(2,d),2)+1;
end
if ~isnumeric(directions)||ndims(directions)>2 ...
        ||size(directions,1)~=d||any(~isfinite(directions(:)))
    error('paraunitary_fir:Directions', ...
        'DIRECTIONS must be finite D-by-M.');
end
directions=double(directions);
for j=1:size(directions,2)
    nv=norm(directions(:,j));
    if nv==0,error('paraunitary_fir:Directions', ...
            'Directions must be nonzero.');end
    directions(:,j)=directions(:,j)/nv;
end
if nargin<4||isempty(U0),U0=eye(d);end
if ~isnumeric(U0)||~isequal(size(U0),[d d]) ...
        ||any(~isfinite(U0(:))) ...
        ||norm(U0'*U0-eye(d),'fro')>5e-11*max(1,d)
    error('paraunitary_fir:U0','U0 must be a finite D-by-D unitary.');
end
if nargin<5||isempty(omega),omega=0;end
validateattributes(omega,{'numeric'},{'scalar','real','finite'});
if nargin<6||isempty(output),output='matrix';end
if ~(ischar(output)||(isstring(output)&&isscalar(output)))
    error('paraunitary_fir:Mode','OUTPUT must be a scalar string.');
end
if ~isfinite(d*L)||d*L>flintmax
    error('paraunitary_fir:Order','D*L exceeds exact index range.');
end
switch lower(char(output))
    case 'frequency'
        U=um_paraunitary_response(exp(1i*omega),directions,U0);
    case 'polyphase'
        U=struct('taps',um_paraunitary_taps(directions,U0), ...
            'winding',-size(directions,2),'channels',d,'period',L, ...
            'delay_exponents',0:size(directions,2));
    case 'spectral'
        U=unitary_spectral_reference('paraunitary_fir',L);
        U.dimension=d*L;
        U.finite.kind='exact';
        U.finite.eigenvalues=@()um_paraunitary_spectral(L,directions,U0,'eigenvalues');
        U.finite.moments=@(k)um_paraunitary_spectral(L,directions,U0,'moments',k);
        U.reference='Vaidyanathan, Multirate Systems and Filter Banks (1993)';
        U.notes=['Exact finite periodic trace spectrum from the L small ' ...
            'D-by-D Fourier blocks, not the spectrum of one response ' ...
            'or a continuous frequency-phase density.'];
    case {'operator','implicit'}
        U=@(X,flag)um_paraunitary_apply(X,flag,L,directions,U0);
    case {'matrix','explicit'}
        if d*L>4096
            error('paraunitary_fir:Order','Use operator/spectral at order >4096.');
        end
        F=@(X,flag)um_paraunitary_apply(X,flag,L,directions,U0);
        U=F(eye(d*L),'notransp');
    otherwise
        error('paraunitary_fir:Mode', ...
            'Use matrix, operator, spectral, frequency, or polyphase.');
end
end
