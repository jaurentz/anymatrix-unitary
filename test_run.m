function results = test_run
%TEST_RUN Run the MATLAB unit tests, standalone or via Anymatrix.
%   RESULTS = TEST_RUN executes the same tests in both environments.
%   Anymatrix clones external groups into unitary/private/. MATLAB does not
%   make private functions available to tests in a nested tests/ directory,
%   even if the private folder is added to the path. To avoid maintaining a
%   second test suite, mirror the MATLAB files in a temporary ordinary
%   folder when called from inside Anymatrix.
root = fileparts(mfilename('fullpath'));
oldPath = path;
temporaryRoot = '';
[~,folderName] = fileparts(root);
if strcmp(folderName,'private')
    temporaryRoot = tempname;
    mkdir(temporaryRoot);
    [ok,msg] = copyfile(fullfile(root,'*.m'),temporaryRoot);
    assert(ok,'Could not copy group M-files to temporary test folder: %s',msg);
    [ok,msg] = copyfile(fullfile(root,'tests'),fullfile(temporaryRoot,'tests'));
    assert(ok,'Could not copy group tests to temporary test folder: %s',msg);
    [ok,msg] = copyfile(fullfile(root,'metadata'),fullfile(temporaryRoot,'metadata'));
    assert(ok,'Could not copy metadata: %s',msg);
    [ok,msg] = copyfile(fullfile(root,'docs'),fullfile(temporaryRoot,'docs'));
    assert(ok,'Could not copy documentation: %s',msg);
    [ok,msg] = copyfile(fullfile(root,'README.md'),fullfile(temporaryRoot,'README.md'));
    assert(ok,'Could not copy README: %s',msg);
    [ok,msg] = copyfile(fullfile(root,'CONTRIBUTING.md'),fullfile(temporaryRoot,'CONTRIBUTING.md'));
    assert(ok,'Could not copy contribution guide: %s',msg);
    root = temporaryRoot;
end
cleanupHandle = onCleanup(@()restore_and_remove(oldPath,temporaryRoot)); %#ok<NASGU>
addpath(root);
results = runtests(fullfile(root,'tests'),'IncludeSubfolders',false);
assertSuccess(results);
end

function restore_and_remove(oldPath,temporaryRoot)
path(oldPath);
if ~isempty(temporaryRoot) && isfolder(temporaryRoot)
    rmdir(temporaryRoot,'s');
end
end
