function [U,properties]=clifford_circuit(q,circuit,mode)
%CLIFFORD_CIRCUIT Finite Clifford/stabilizer circuits with local-gate actions.
%   U=CLIFFORD_CIRCUIT(Q,CIRCUIT,MODE) returns an exactly specified
%   2^Q-by-2^Q Clifford unitary. Q is the qubit count (1<=Q<=20);
%   qubit 1 is the most-significant bit in |b1...bQ>.
%   CIRCUIT is a preset name, N-by-2 cell {NAME,SITES}, or ordered
%   struct array with fields NAME and SITES.
%   Gates are applied chronologically by rows: U=Gm*...*G2*G1.
%   Default: 'hadamard' for Q=1; 'bell' for Q>=2.
%
%   Presets: identity, hadamard, phase, phase_dagger, pauli_x,
%   pauli_y, pauli_z, cnot, cz, swap, bell.
%   Custom names: H,S,SDG,X,Y,Z,CNOT,CX,CZ,SWAP with ordered sites.
%   There is NO free scalar global phase: H=[1 1;1 -1]/sqrt(2),
%   S=diag(1,i), SDG=diag(1,-i), Y=[0 -i;i 0]; CNOT uses the
%   first listed site as control and the standard local basis.
%
%   MODE='matrix' (default) gives dense 2^Q unitary, Q<=10.
%   MODE='operator' gives O(X,FLAG) using the shared local-gate engine.
%   MODE='spectral' gives exact finite trace phase multiplicities
%   for identity and single-gate circuits; generic words' analytic
%   finite spectra are marked UNKNOWN. These gate words have finite
%   order under the fixed phase convention, but NOT a CUE law.
%
%   J. Dehaene and B. De Moor, Phys. Rev. A 68 (2003), 042318 [128].
%     https://doi.org/10.1103/PhysRevA.68.042318
%   D. Gross, S. Nezami and M. Walter, Commun. Math. Phys. 385
%   (2021), 1325-1393 [129].
%     https://doi.org/10.1007/s00220-021-04118-7
properties={'square','unitary','normal','complex', ...
    'parameter-dependent','scalable'};
if nargin==0,U=[];return;end
validateattributes(q,{'numeric'},{'scalar','real','finite','integer','>=',1,'<=',20});
if nargin<2 || isempty(circuit)
    if q==1,circuit='hadamard';else,circuit='bell';end
end
if nargin<3 || isempty(mode),mode='matrix';end
if ~(ischar(mode)||(isstring(mode)&&isscalar(mode)))
    error('clifford_circuit:Mode','MODE must be text.');
end
mode=lower(char(mode));
if ~any(strcmp(mode,{'matrix','operator','spectral'}))
    error('clifford_circuit:Mode','Unknown output mode.');
end
list=parse_circuit(q,circuit);
ng=size(list,1);
gates=struct('gate',{},'sites',{});
names=cell(ng,1);
for t=1:ng
    [name,gate,sites]=compile_gate(q,list{t,1},list{t,2});
    gates(t)=struct('gate',gate,'sites',sites); %#ok<AGROW>
    names{t}=name;
end
if ng==0
    profile='identity';
elseif ng==1
    profile=lower(names{1});
else
    profile='unknown';
end
if strcmp(mode,'spectral')
    U=um_clifford_spectral(q,profile);
    return
end
if strcmp(mode,'matrix') && q>10
    error('clifford_circuit:ExplicitOrder', ...
        'Dense output supports Q<=10; use operator or spectral mode.');
end
N=2^q;
O=@(X,flag)um_circuit_apply(X,flag,q,gates,1);
if strcmp(mode,'operator')
    U=O;
else
    U=O(eye(N),'notransp');
end
end

function list=parse_circuit(q,circuit)
if ischar(circuit)||(isstring(circuit)&&isscalar(circuit))
    name=lower(strtrim(char(circuit)));
    switch name
        case 'identity'
            list=cell(0,2);
        case 'hadamard'
            list={'H',1};
        case 'phase'
            list={'S',1};
        case 'phase_dagger'
            list={'SDG',1};
        case 'pauli_x'
            list={'X',1};
        case 'pauli_y'
            list={'Y',1};
        case 'pauli_z'
            list={'Z',1};
        case 'cnot'
            list={'CNOT',[1 2]};
        case 'cz'
            list={'CZ',[1 2]};
        case 'swap'
            list={'SWAP',[1 2]};
        case 'bell'
            list={'H',1;'CNOT',[1 2]};
        otherwise
            error('clifford_circuit:Circuit','Unknown Clifford preset.');
    end
elseif iscell(circuit)
    if isempty(circuit)
        list=cell(0,2);
    elseif ndims(circuit)==2 && size(circuit,2)==2
        list=circuit;
    else
        error('clifford_circuit:Circuit', ...
            'CIRCUIT must be an N-by-2 cell array of {name,sites}.');
    end
elseif isstruct(circuit)
    if ~all(isfield(circuit,{'name','sites'}))
        error('clifford_circuit:Circuit','Struct gates require name and sites.');
    end
    list=cell(numel(circuit),2);
    for k=1:numel(circuit)
        list{k,1}=circuit(k).name;
        list{k,2}=circuit(k).sites;
    end
else
    error('clifford_circuit:Circuit','Invalid circuit description.');
end
if q==1 && ~isempty(list)
    for j=1:size(list,1)
        if (ischar(list{j,1}) || isstring(list{j,1})) && ...
                any(strcmpi(char(list{j,1}),{'CNOT','CX','CZ','SWAP'}))
            error('clifford_circuit:Sites', ...
                'Two-qubit gates require at least two qubits.');
        end
    end
end
end

function [name,gate,sites]=compile_gate(q,rawName,sites)
if ~(ischar(rawName)||(isstring(rawName)&&isscalar(rawName)))
    error('clifford_circuit:Gate','Gate name must be text.');
end
name=upper(strtrim(char(rawName)));
if strcmp(name,'CX'),name='CNOT';end
if ~(isnumeric(sites) && isreal(sites) && isvector(sites) ...
        && ~isempty(sites) && all(isfinite(sites)) ...
        && all(sites==round(sites)) && all(sites>=1) ...
        && all(sites<=q) && numel(unique(sites))==numel(sites))
    error('clifford_circuit:Sites','Sites must be distinct valid qubits.');
end
sites=double(sites(:).');
switch name
    case 'H'
        gate=[1 1;1 -1]/sqrt(2);
    case 'S'
        gate=diag([1,1i]);
    case 'SDG'
        gate=diag([1,-1i]);
    case 'X'
        gate=[0 1;1 0];
    case 'Y'
        gate=[0 -1i;1i 0];
    case 'Z'
        gate=diag([1,-1]);
    case 'CNOT'
        gate=[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];
    case 'CZ'
        gate=diag([1 1 1 -1]);
    case 'SWAP'
        gate=[1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1];
    otherwise
        error('clifford_circuit:Gate','Unknown Clifford gate.');
end
if numel(sites)~=round(log2(size(gate,1)))
    error('clifford_circuit:Sites','Gate arity does not match sites.');
end
end
