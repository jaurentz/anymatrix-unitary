function tests=test_catalog
%TEST_CATALOG Public inventory, documentation, and property-list consistency.
tests=functiontests(localfunctions);
end

function test_inventory_matches_generators(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
raw=fileread(fullfile(root,'metadata','catalog.json'));
meta=jsondecode(raw);
entries=meta.matrix_families;
verifyEqual(testCase,numel(entries),22);
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
end
end

function test_public_documents_are_self_contained(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
files={'README.md','docs/MATHEMATICS.md','docs/MATRIX_ENCYCLOPEDIA.md', ...
       'docs/TESTING.md','CONTRIBUTING.md'};
for k=1:numel(files)
    document=fileread(fullfile(root,files{k}));
    verifyGreaterThan(testCase,numel(document),200);
end
encyclopedia=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyTrue(testCase,contains(encyclopedia,'## Bibliography'));
verifyTrue(testCase,contains(encyclopedia,'## General conventions'));
maths=fileread(fullfile(root,'docs','MATHEMATICS.md'));
verifyFalse(testCase,contains(maths,'\['));
verifyFalse(testCase,contains(maths,'\]'));
verifyGreaterThan(testCase,count(maths,'$$'),6);
end
end
