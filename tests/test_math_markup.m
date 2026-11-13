function tests=test_math_markup
%TEST_MATH_MARKUP Check Markdown mathematics conventions in source.
tests=functiontests(localfunctions);
end

function test_display_math_fences(tc)
root=fileparts(fileparts(mfilename('fullpath')));
paths={'docs/MATRIX_ENCYCLOPEDIA.md','docs/MATHEMATICS.md','README.md'};
minimum=[15 10 1];
for j=1:numel(paths)
    s=fileread(fullfile(root,paths{j}));
    starts=regexp(s,'(?m)^```math\r?$','match');
    blocks=regexp(s,'(?ms)^```math\r?\n(.*?)\r?\n```\r?$','tokens');
    verifyGreaterThanOrEqual(tc,numel(starts),minimum(j),paths{j});
    verifyEqual(tc,numel(blocks),numel(starts),paths{j});
    verifyFalse(tc,contains(s,repmat(char(36),1,2)),paths{j});
    verifyFalse(tc,any(strtrim(splitlines(string(s)))=="$"),paths{j});
    for k=1:numel(blocks)
        verifyNotEmpty(tc,strtrim(blocks{k}{1}),paths{j});
    end
end
end

function test_cazac_display_math(tc)
root=fileparts(fileparts(mfilename('fullpath')));
s=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
a=strfind(s,'### Chu/CAZAC unitary circulants');
b=strfind(s,'### Discrete Fresnel transform');
verifyEqual(tc,numel(a),1);verifyEqual(tc,numel(b),1);
section=s(a:b-1);
verifyEqual(tc,numel(regexp(section,'(?m)^```math\r?$','match')),3);
verifyFalse(tc,contains(section,'\begin{cases}'));
verifyFalse(tc,contains(section,'\bmod'));
verifyTrue(tc,contains(section,'c_k=\exp'));
verifyTrue(tc,contains(section,'U_{jk}=\frac{c_d}'));
end

function test_julia_halmos_display_math(tc)
root=fileparts(fileparts(mfilename('fullpath')));
s=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
a=strfind(s,'### 22. Julia/Halmos unitary extension');
b=strfind(s,'## F. Further real and complex OPUC');
verifyEqual(tc,numel(a),1);verifyEqual(tc,numel(b),1);
section=s(a:b-1);
verifyEqual(tc,numel(regexp(section,'(?m)^```math\r?$','match')),4);
verifyTrue(tc,contains(section,'\mathcal{U}(T)'));
verifyTrue(tc,contains(section,'D_{T^*}'));
verifyTrue(tc,contains(section,'P\mathcal{U}(T)^kP^*'));
end

function test_indexing_and_terminology(tc)
root=fileparts(fileparts(mfilename('fullpath')));
s=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyTrue(tc,contains(s,'e_{n+1}=e_1'));
verifyTrue(tc,contains(s,'\mathbf{1}_{n\mid k}'));
verifyTrue(tc,contains(s,'\mathcal{U}(T)'));
readme=fileread(fullfile(root,'README.md'));
verifyTrue(tc,contains(readme,'\alpha_{2k}'));
verifyTrue(tc,contains(readme,'\alpha_{2k+1}'));
verifyTrue(tc,contains(readme,'examples/unit_circle_spectral_bands.m'));
end
