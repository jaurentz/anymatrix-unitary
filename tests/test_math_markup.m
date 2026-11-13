function tests = test_math_markup
%TEST_MATH_MARKUP Check GitHub-compatible mathematical Markdown conventions.
%   This is a source-format regression check; it does not replace visual
%   inspection of GitHub's rendered mathematics.
tests = functiontests(localfunctions);
end

function test_display_math_fences(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
paths = {'docs/MATRIX_ENCYCLOPEDIA.md','docs/MATHEMATICS.md'};
minimum = [15,10];
for j = 1:numel(paths)
    source = fileread(fullfile(root,paths{j}));
    starts = regexp(source,'(?m)^```math\r?$','match');
    blocks = regexp(source,'(?ms)^```math\r?\n(.*?)\r?\n```\r?$','tokens');
    verifyGreaterThanOrEqual(testCase,numel(starts),minimum(j),paths{j});
    verifyEqual(testCase,numel(blocks),numel(starts),paths{j});
    verifyFalse(testCase,contains(source,repmat(char(36),1,2)),paths{j});
    verifyFalse(testCase,any(strtrim(splitlines(string(source)))=="$"),paths{j});
    verifyFalse(testCase,contains(source,'\operatorname'),paths{j});
    verifyFalse(testCase,contains(source,'\DeclareMathOperator'),paths{j});
    verifyFalse(testCase,contains(source,'\pmod'),paths{j});
    verifyFalse(testCase,contains(source,'\mathbf1'),paths{j});
    verifyFalse(testCase,contains(source,'\overline\alpha'),paths{j});
    verifyFalse(testCase,contains(source,'\['),paths{j});
    verifyFalse(testCase,contains(source,'\]'),paths{j});
    for k = 1:numel(blocks)
        verifyNotEmpty(testCase,strtrim(blocks{k}{1}),paths{j});
    end
end
end

function test_indexing_and_terminology(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyTrue(testCase,contains(source,'e_{n+1}=e_1'));
verifyTrue(testCase,contains(source,'\mathbf{1}_{n\mid k}'));
verifyTrue(testCase,contains(source,'\mathcal{U}(T)'));
end
