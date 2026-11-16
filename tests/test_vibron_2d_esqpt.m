function tests=test_vibron_2d_esqpt
%TEST_VIBRON_2D_ESQPT Independent U(3) boson reference and spectral tests.
tests=functiontests(localfunctions);
end

function test_boson_operators_independent_construction(tc)
% Independently assemble Dplus from sigma and tau+/- creation and
% annihilation operators, not from the tridiagonal sector formula.
for N=[1 2 3 4 6 8]
    d=(N+1)*(N+2)/2;
    states=zeros(d,3);ells=zeros(d,1);ts=zeros(d,1);
    index=zeros(N+1,N+1);
    pos=0;
    for ell=-N:N
        for tau=abs(ell):2:N
            pos=pos+1;
            p=(tau+ell)/2;m=(tau-ell)/2;s=N-tau;
            states(pos,:)=[p,m,s];
            ells(pos)=ell;ts(pos)=tau;
            index(p+1,m+1)=pos;
        end
    end
    verifyEqual(tc,pos,d);
    rows=[];cols=[];vals=[];
    for j=1:d
        p=states(j,1);m=states(j,2);s=states(j,3);
        if s>0
            rows(end+1)=index(p+2,m+1); %#ok<AGROW>
            cols(end+1)=j; %#ok<AGROW>
            vals(end+1)=sqrt(2*s*(p+1)); %#ok<AGROW>
        end
        if m>0
            rows(end+1)=index(p+1,m); %#ok<AGROW>
            cols(end+1)=j; %#ok<AGROW>
            vals(end+1)=-sqrt(2*m*(s+1)); %#ok<AGROW>
        end
    end
    Dplus=sparse(rows,cols,vals,d,d);
    W2=(Dplus*Dplus'+Dplus'*Dplus)/2+ ...
        spdiags(ells.^2,0,d,d);
    for xi=[0 .13 .5 1]
        href=(1-xi)*spdiags(ts,0,d,d)+ ...
            (xi/N)*(N*(N+1)*speye(d)-W2);
        H=vibron_2d_esqpt(N,xi,.22,[],'hamiltonian');
        verifyTrue(tc,issparse(H));
        verifyEqual(tc,size(H),[d d]);
        verifyLessThan(tc,norm(full(H-href),'fro'),5e-12*d);
        verifyLessThan(tc,norm(full(H-H'),'fro'),3e-12*d);
    end
end
end

function test_exact_dynamical_symmetries_and_moments(tc)
for N=[1 2 4 7 10]
    for ell=num2cell([NaN 0 1 -1 N])
        l=ell{1};if isnan(l),l=[];end
        for xi=[0 1]
            H=vibron_2d_esqpt(N,xi,.41,l,'hamiltonian');
            S=vibron_2d_esqpt(N,xi,.41,l,'spectral');
            U=vibron_2d_esqpt(N,xi,.41,l);
            z=S.finite.eigenvalues();
            verifyEqual(tc,S.finite.kind,'exact');
            verifyEqual(tc,numel(z),size(H,1));
            verifyLessThan(tc,norm(full(U'*U-speye(size(U,1))),'fro'), ...
                5e-12*size(U,1));
            e=eig(full(H));
            actual=exp(-1i*.41*e);
            verifyLessThan(tc,hausdorff(actual,z), ...
                5e-12*max(1,size(H,1)));
            for k=[-2 0 1 3]
                verifyLessThan(tc,abs(trace(U^k)/numel(z)-S.finite.moments(k)), ...
                    1e-11*max(1,numel(z)));
            end
        end
    end
end
end

function test_direct_sum_and_plus_minus_angular_momentum(tc)
N=18;xi=.42;time=.13;
H=vibron_2d_esqpt(N,xi,time,[],'hamiltonian');
U=vibron_2d_esqpt(N,xi,time);
d=(N+1)*(N+2)/2;
verifyEqual(tc,size(U),[d d]);
verifyLessThan(tc,norm(full(U'*U-speye(d)),'fro'),1e-10);
blocks=cell(2*N+1,1);Hblocks=blocks;
for ell=-N:N
    idx=ell+N+1;
    blocks{idx}=vibron_2d_esqpt(N,xi,time,ell);
    Hblocks{idx}=vibron_2d_esqpt(N,xi,time,ell,'hamiltonian');
    verifyLessThan(tc,norm(full(Hblocks{idx}- ...
       vibron_2d_esqpt(N,xi,time,-ell,'hamiltonian')),'fro'), ...
       5e-12);
end
verifyLessThan(tc,norm(full(U-blkdiag(blocks{:})),'fro'),1e-10);
verifyLessThan(tc,norm(full(H-blkdiag(Hblocks{:})),'fro'),1e-10);
end

function test_algebraically_shrinking_sector_gap(tc)
% Fig. 2(b) of Khalouf-Rivera et al. (2024), low-lying ell=0,1 pairs.
sizes=[20 40 80];gaps=zeros(size(sizes));
for j=1:numel(sizes)
    N=sizes(j);
    E0=eig(full(vibron_2d_esqpt(N,.5,1,0,'hamiltonian')));
    E1=eig(full(vibron_2d_esqpt(N,.5,1,1,'hamiltonian')));
    gaps(j)=abs(min(real(E0))-min(real(E1)));
end
verifyTrue(tc,all(diff(gaps)<0));
verifyGreaterThan(tc,gaps(1)/gaps(2),1.5);
verifyGreaterThan(tc,gaps(2)/gaps(3),1.5);
verifyGreaterThan(tc,gaps(3),1e-4);
end

function test_unknown_spectrum_modes_and_errors(tc)
S=vibron_2d_esqpt(30,.5,.16,1,'spectral');
verifyEqual(tc,S.family,'vibron_2d_esqpt');
verifyEqual(tc,S.finite.kind,'unknown');
verifyEqual(tc,S.dimension,15);
T=vibron_2d_esqpt(12,.5,0,[],'spectral');
verifyEqual(tc,T.finite.kind,'exact');
verifyEqual(tc,T.finite.eigenvalues(),ones(91,1));
[empty,properties]=vibron_2d_esqpt();
verifyEmpty(tc,empty);
verifyTrue(tc,any(strcmp(properties,'unitary')));
verifyError(tc,@()vibron_2d_esqpt(6,-.1), ...
    'MATLAB:notGreaterEqual');
verifyError(tc,@()vibron_2d_esqpt(6,.5,.1,7), ...
    'MATLAB:notLessEqual');
verifyError(tc,@()vibron_2d_esqpt(6,.5,.1,[],'bogus'), ...
    'vibron_2d_esqpt:Mode');
verifyError(tc,@()vibron_2d_esqpt(129,.5,.1), ...
    'vibron_2d_esqpt:ExplicitOrder');
end

function d=hausdorff(a,b)
a=a(:);b=b(:);
d=max([max(arrayfun(@(v)min(abs(v-b)),a)), ...
       max(arrayfun(@(v)min(abs(v-a)),b))]);
end
