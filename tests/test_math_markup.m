function tests = test_math_markup
%TEST_MATH_MARKUP Check GitHub-compatible mathematical Markdown conventions.
%   This is a source-format regression check; it does not replace visual
%   inspection of GitHub's rendered mathematics.
tests = functiontests(localfunctions);
end

function test_display_math_fences(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
paths = {'docs/MATRIX_ENCYCLOPEDIA.md','docs/MATHEMATICS.md','README.md'};
minimum = [15,10,1];
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
    verifyFalse(testCase,contains(source,'alpha_(2k)'),paths{j});
    verifyFalse(testCase,contains(source,'z^(-1)'),paths{j});
    for k = 1:numel(blocks)
        verifyNotEmpty(testCase,strtrim(blocks{k}{1}),paths{j});
    end
end
end

function test_cazac_display_math(testCase)
% The two Zadoff-Chu parity formulas and circulant definition are kept
% in separate GitHub math fences rather than a combined cases expression.
root = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
startIndex = strfind(source,'### Chu/CAZAC unitary circulants');
endIndex = strfind(source,'### Discrete Fresnel transform');
verifyEqual(testCase,numel(startIndex),1);
verifyEqual(testCase,numel(endIndex),1);
section = source(startIndex:endIndex-1);
verifyEqual(testCase,numel(regexp(section,'(?m)^```math
% The Julia/Halmos section uses fenced mathematics to avoid raw dollar signs
% appearing in the rendered GitHub encyclopedia.
root = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
startIndex = strfind(source,'### 22. Julia/Halmos unitary extension');
endIndex = strfind(source,'## F. Further real and complex OPUC');
verifyEqual(testCase,numel(startIndex),1);
verifyEqual(testCase,numel(endIndex),1);
section = source(startIndex:endIndex-1);
verifyFalse(testCase,contains(section,char(36)));
verifyEqual(testCase,numel(regexp(section,'(?m)^```math
root = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyTrue(testCase,contains(source,'e_{n+1}=e_1'));
verifyTrue(testCase,contains(source,'\mathbf{1}_{n\mid k}'));
verifyTrue(testCase,contains(source,'\mathcal{U}(T)'));
readme = fileread(fullfile(root,'README.md'));
verifyTrue(testCase,contains(readme,'\alpha_{2k}'));
verifyTrue(testCase,contains(readme,'\alpha_{2k+1}'));
verifyTrue(testCase,contains(readme,'examples/unit_circle_spectral_bands.m'));
verifyFalse(testCase,contains(readme,'showcase_new_families'));
end
,'match')),4);
verifyTrue(testCase,contains(section,'\mathcal{U}(T)'));
verifyTrue(testCase,contains(section,'D_{T^*}'));
verifyTrue(testCase,contains(section,'P\mathcal{U}(T)^kP^*'));
end

function test_indexing_and_terminology(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyTrue(testCase,contains(source,'e_{n+1}=e_1'));
verifyTrue(testCase,contains(source,'\mathbf{1}_{n\mid k}'));
verifyTrue(testCase,contains(source,'\mathcal{U}(T)'));
readme = fileread(fullfile(root,'README.md'));
verifyTrue(testCase,contains(readme,'\alpha_{2k}'));
verifyTrue(testCase,contains(readme,'\alpha_{2k+1}'));
verifyTrue(testCase,contains(readme,'examples/unit_circle_spectral_bands.m'));
verifyFalse(testCase,contains(readme,'showcase_new_families'));
end
,'match')),3);
verifyFalse(testCase,contains(section,'\begin{cases}'));
verifyFalse(testCase,contains(section,'\bmod'));
verifyTrue(testCase,contains(section,'c_k=\exp'));
verifyTrue(testCase,contains(section,'U_{jk}=\frac{c_d}'));
end

function test_julia_halmos_display_math(testCase)
% The Julia/Halmos section uses fenced mathematics to avoid raw dollar signs
% appearing in the rendered GitHub encyclopedia.
root = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
startIndex = strfind(source,'### 22. Julia/Halmos unitary extension');
endIndex = strfind(source,'## F. Further real and complex OPUC');
verifyEqual(testCase,numel(startIndex),1);
verifyEqual(testCase,numel(endIndex),1);
section = source(startIndex:endIndex-1);
verifyFalse(testCase,contains(section,char(36)));
verifyEqual(testCase,numel(regexp(section,'(?m)^```math
root = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyTrue(testCase,contains(source,'e_{n+1}=e_1'));
verifyTrue(testCase,contains(source,'\mathbf{1}_{n\mid k}'));
verifyTrue(testCase,contains(source,'\mathcal{U}(T)'));
readme = fileread(fullfile(root,'README.md'));
verifyTrue(testCase,contains(readme,'\alpha_{2k}'));
verifyTrue(testCase,contains(readme,'\alpha_{2k+1}'));
verifyTrue(testCase,contains(readme,'examples/unit_circle_spectral_bands.m'));
verifyFalse(testCase,contains(readme,'showcase_new_families'));
end
,'match')),4);
verifyTrue(testCase,contains(section,'\mathcal{U}(T)'));
verifyTrue(testCase,contains(section,'D_{T^*}'));
verifyTrue(testCase,contains(section,'P\mathcal{U}(T)^kP^*'));
end

function test_indexing_and_terminology(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
source = fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyTrue(testCase,contains(source,'e_{n+1}=e_1'));
verifyTrue(testCase,contains(source,'\mathbf{1}_{n\mid k}'));
verifyTrue(testCase,contains(source,'\mathcal{U}(T)'));
readme = fileread(fullfile(root,'README.md'));
verifyTrue(testCase,contains(readme,'\alpha_{2k}'));
verifyTrue(testCase,contains(readme,'\alpha_{2k+1}'));
verifyTrue(testCase,contains(readme,'examples/unit_circle_spectral_bands.m'));
verifyFalse(testCase,contains(readme,'showcase_new_families'));
end
