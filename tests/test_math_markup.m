function tests=test_math_markup
%TEST_MATH_MARKUP Check Markdown mathematics conventions in source.
tests=functiontests(localfunctions);
end

function test_display_math_fences(tc)
root=fileparts(fileparts(mfilename('fullpath')));
paths={'docs/MATRIX_ENCYCLOPEDIA.md','docs/MATHEMATICS.md','README.md', ...
       'docs/STRUCTURED_OPERATORS.md','docs/SPECTRAL_REFERENCE.md', ...
       'docs/SPECTRAL_DISTRIBUTION_ATLAS.md'};
minimum=[15 10 1 9 1 2];
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

function test_affine_fourier_equation_renders(tc)
% Componentwise chirp notation should stay compatible with GitHub math.
root=fileparts(fileparts(mfilename('fullpath')));
s=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
a=strfind(s,'### Affine Fourier / AFDM modulation');
b=strfind(s,'### Finite Zak transforms');
verifyEqual(tc,numel(a),1);
verifyEqual(tc,numel(b),1);
section=s(a:b-1);
verifyEqual(tc,numel(strfind(section,'```math')),1);
verifyTrue(tc,contains(section,'D(c)_{kk}=e^{-2'));
verifyFalse(tc,contains(section,'\mathop'));
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

function test_all_public_math_notation(tc)
% Check GitHub-compatible delimiters, macros and single-line inline math.
root=fileparts(fileparts(mfilename('fullpath')));
paths={'README.md','docs/MATHEMATICS.md','docs/MATRIX_ENCYCLOPEDIA.md', ...
       'docs/STRUCTURED_OPERATORS.md','docs/SPECTRAL_REFERENCE.md', ...
       'docs/SPECTRAL_DISTRIBUTION_ATLAS.md'};
for j=1:numel(paths)
    s=fileread(fullfile(root,paths{j}));
    verifyFalse(tc,contains(s,'\operatorname'),paths{j});
    verifyFalse(tc,contains(s,'\DeclareMathOperator'),paths{j});
    verifyFalse(tc,contains(s,'\newcommand'),paths{j});
    verifyFalse(tc,contains(s,'\tag{'),paths{j});
    verifyEmpty(tc,regexp(s,'(?<![A-Za-z\\])o\s*\(\s*[Nn]','match'), ...
        sprintf('Avoid little-o complexity notation in %s',paths{j}));
    verifyFalse(tc,contains(s,'\('),paths{j});
    verifyFalse(tc,contains(s,'\)'),paths{j});
    verifyEmpty(tc,regexp(s,'(?m)^\s*\\\[\s*$','match'),paths{j});
    verifyEmpty(tc,regexp(s,'(?m)^\s*\\\]\s*$','match'),paths{j});
    lines=splitlines(string(s));
    inFence=false;
    for k=1:numel(lines)
        line=strtrim(lines(k));
        if startsWith(line,"```") || startsWith(line,"~~~")
            inFence=~inFence;
            continue
        end
        if inFence
            continue
        end
        dollars=regexp(char(lines(k)),'(?<!\\)\$','match');
        verifyEqual(tc,mod(numel(dollars),2),0, ...
            sprintf('Unpaired inline math delimiter in %s line %d',paths{j},k));
        % Verify that inline math is not directly glued to letter suffixes.
        % GitHub Markdown can display $n$th literally rather than render n.
        lineText=char(lines(k));
        mathEnds=regexp(lineText,'(?<!\\)\$[^$\n]+\$','end');
        for m=1:numel(mathEnds)
            if mathEnds(m)<numel(lineText)
                verifyFalse(tc,isletter(lineText(mathEnds(m)+1)), ...
                    sprintf('Math suffix adjacent to text in %s line %d',paths{j},k));
            end
        end
        % Avoid fragile hyphenated formula headings, e.g. circular-$\beta$.
        verifyEmpty(tc,regexp(lineText,'[A-Za-z]-\$\\beta\$','match'), ...
            sprintf('Hyphenated beta math in %s line %d',paths{j},k));
    end
end
end

function test_matrix_order_convention(tc)
% All explicit encyclopedia order declarations use the actual N-by-N order.
root=fileparts(fileparts(mfilename('fullpath')));
ency=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
orders=regexp(ency,'(?m)^\*\*Order:\*\*[^\n]*','match');
verifyGreaterThan(tc,numel(orders),80);
for k=1:numel(orders)
    verifyTrue(tc,startsWith(orders{k},'**Order:** $N'), ...
        sprintf('Matrix-order notation in encyclopedia declaration %d',k));
end
verifyTrue(tc,contains(ency,'U\in\mathbb C^{N\times N}'));
structure=fileread(fullfile(root,'docs','STRUCTURED_OPERATORS.md'));
verifyTrue(tc,contains(structure,'$N\times N$'));
matrixFree=fileread(fullfile(root,'docs','MATRIX_FREE.md'));
verifyTrue(tc,contains(matrixFree,'N=2*n'));
end
