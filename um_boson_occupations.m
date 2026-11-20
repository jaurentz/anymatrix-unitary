function occ=um_boson_occupations(d,p)
%UM_BOSON_OCCUPATIONS Deterministic normalized Fock basis occupations.
%   The ascending separators enumerate weak D-part compositions of P.
if d==1,occ=p;return;end
separators=nchoosek(1:(d+p-1),d-1);
N=size(separators,1);
occ=diff([zeros(N,1),separators,(d+p)*ones(N,1)],1,2)-1;
end
