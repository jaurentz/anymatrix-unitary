function tests=test_catalog
%TEST_CATALOG Public inventory, documentation, and property-list consistency.
tests=functiontests(localfunctions);
end

function test_inventory_matches_generators(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
raw=fileread(fullfile(root,'metadata','catalog.json'));
meta=jsondecode(raw);
entries=meta.matrix_families;
verifyEqual(testCase,numel(entries),138);
ids={entries.id};
verifyEqual(testCase,numel(unique(ids)),numel(ids));
encyclopedia=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
recognized={'banded','bidiagonal','binary','block Toeplitz','built-in', ...
'circulant','complex','correlation','defective','diagonally dominant', ...
'eigenvalues','fixed size','hankel','hermitian','hessenberg','idempotent', ...
'indefinite','ill conditioned','infinitely divisible','integer','inverse', ...
'involutory','M-matrix','nilpotent','nonnegative','normal','orthogonal', ...
'parameter-dependent','permutation','positive','positive definite', ...
'pseudo-orthogonal','random','rank deficient','real','rectangular', ...
'scalable','singular values','skew-hermitian','skew-symmetric','sparse', ...
'square','stochastic','symmetric','toeplitz','totally nonnegative', ...
'totally positive','triangular','tridiagonal','unimodular','unitary'};
for k=1:numel(entries)
    name=entries(k).id;
    verifyTrue(testCase,isfile(fullfile(root,[name,'.m'])), ...
        sprintf('Missing generator %s',name));
    verifyTrue(testCase,contains(encyclopedia,['`',name,'.m`']), ...
        sprintf('Missing encyclopedia entry for %s',name));
    [A,props]=feval(name);
    verifyEmpty(testCase,A);
    verifyTrue(testCase,iscellstr(props)); %#ok<ISCLSTR>
    verifyTrue(testCase,all(ismember(props,recognized)), ...
        sprintf('Unrecognized Anymatrix property for %s',name));
    verifyTrue(testCase,any(strcmp(props,'unitary')));
    entry=entries(k);
    verifyTrue(testCase,islogical(entry.matrix_free_supported));
    if entry.matrix_free_supported
        kernels=strsplit(char(entry.operator_kernel),' / ');
        for ik=1:numel(kernels)
            kernel=strtrim(kernels{ik});
            verifyTrue(testCase,startsWith(kernel,'um_'), ...
                sprintf('Invalid matrix-free kernel %s for %s',kernel,name));
            verifyTrue(testCase,isfile(fullfile(root,[kernel,'.m'])), ...
                sprintf('Missing matrix-free kernel %s for %s',kernel,name));
        end
    end
end
end

function test_public_documents_are_self_contained(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
files={'README.md','docs/MATHEMATICS.md','docs/MATRIX_ENCYCLOPEDIA.md', ...
       'docs/TESTING.md','CONTRIBUTING.md'};
for k=1:numel(files)
    document=fileread(fullfile(root,files{k}));
    verifyGreaterThan(testCase,numel(document),200);
    if any(strcmp(files{k},{'docs/MATHEMATICS.md','docs/MATRIX_ENCYCLOPEDIA.md'}))
        verifyFalse(testCase,contains(document,'\operatorname'), ...
            sprintf('Unsupported GitHub math macro in %s',files{k}));
        verifyFalse(testCase,contains(document,'\DeclareMathOperator'), ...
            sprintf('Unsupported GitHub math declaration in %s',files{k}));
    end
end
encyclopedia=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyTrue(testCase,contains(encyclopedia,'## Bibliography'));
verifyTrue(testCase,contains(encyclopedia,'## General conventions'));
verifyTrue(testCase,contains(encyclopedia,'```math'));
verifyTrue(testCase,contains(encyclopedia,'\lVert T\rVert'));
verifyFalse(testCase,contains(encyclopedia,'\mathcal U(T)'));
maths=fileread(fullfile(root,'docs','MATHEMATICS.md'));
verifyFalse(testCase,contains(maths,'\['));
verifyFalse(testCase,contains(maths,'\]'));
verifyGreaterThan(testCase,count(maths,'```math'),6);
end

function test_dirac_catalog_operator_contract(tc)
% Guard against documenting an operator that the catalog marks absent.
root=fileparts(fileparts(mfilename('fullpath')));
meta=jsondecode(fileread(fullfile(root,'metadata','catalog.json')));
names={meta.matrix_families.id};
idx=find(strcmp(names,'dirac_automaton'));
verifyEqual(tc,numel(idx),1);
entry=meta.matrix_families(idx);
verifyTrue(tc,entry.matrix_free_supported);
verifyEqual(tc,char(entry.operator_kernel),'um_dirac_apply');
A=dirac_automaton(7,.35);
U=dirac_automaton(7,.35,'operator');
X=complex(reshape(1:42,14,3),reshape(42:-1:1,14,3))/50;
verifyEqual(tc,U([],'size'),[14 14]);
verifyLessThan(tc,norm(U(X,'notransp')-A*X,'fro'),1e-12);
verifyLessThan(tc,norm(U(X,'transp')-A'*X,'fro'),1e-12);
verifyLessThan(tc,norm(U(X,'ctransp')-A'*X,'fro'),1e-12);
end

function test_documented_inventory_does_not_freeze_stale_counts(tc)
root=fileparts(fileparts(mfilename('fullpath')));
readme=fileread(fullfile(root,'README.md'));
verifyTrue(tc,contains(readme,'metadata/catalog.json'));
verifyEmpty(tc,regexp(readme,'\b[0-9]+\s+(public\s+)?generators\b','match'));
end
