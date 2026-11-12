function results = test_run
%TEST_RUN Run all local unit tests (also called by Anymatrix test command).
%   TEST_RUN is the Anymatrix group test hook. It discovers the MATLAB
%   function-based tests in tests/ and fails if any test fails.
%   RESULTS = TEST_RUN returns matlab.unittest.TestResult values.
root = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(root);
results = runtests(fullfile(root,'tests'),'IncludeSubfolders',false);
assertSuccess(results);
end
