function tests=test_reference_registry
%TEST_REFERENCE_REGISTRY Offline integrity of bibliography and citations.
%   Live DOI resolution was checked outside CI to avoid network-dependent tests.
tests=functiontests(localfunctions);
end

function test_all_bibliography_entries_have_identifiers(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
catalog=jsondecode(fileread(fullfile(root,'metadata','catalog.json')));
registry=jsondecode(fileread(fullfile(root,'metadata','references.json')));
encyclopedia=fileread(fullfile(root,'docs','MATRIX_ENCYCLOPEDIA.md'));
verifyEqual(testCase,catalog.reference_registry,'metadata/references.json');
verifyEqual(testCase,numel(registry.references),46);
for number=1:46
    ref=registry.references(number);
    verifyEqual(testCase,ref.number,number);
    verifyTrue(testCase,contains(encyclopedia,sprintf('%d. ',number)));
    verifyTrue(testCase,startsWith(ref.permalink,'https://'));
    dois=as_cell(ref.dois);
    if isempty(dois)
        verifyTrue(testCase,contains(encyclopedia,ref.permalink));
    else
        for j=1:numel(dois)
            doi=dois{j};
            verifyFalse(testCase,isempty(regexp(doi,'^10\\.[0-9]{4,9}/\\S+$','once')));
            url=['https://doi.org/',doi];
            verifyTrue(testCase,contains(lower(encyclopedia),lower(url)), ...
                sprintf('Unlisted DOI in bibliography reference %d',number));
        end
    end
end
end

function test_all_generator_reference_identifiers(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
catalog=jsondecode(fileread(fullfile(root,'metadata','catalog.json')));
registry=jsondecode(fileread(fullfile(root,'metadata','references.json')));
for k=1:numel(catalog.matrix_families)
    e=catalog.matrix_families(k);
    source=lower(fileread(fullfile(root,[e.id,'.m'])));
    ids=e.encyclopedia_references;
    verifyNotEmpty(testCase,ids,sprintf('No bibliography IDs in %s',e.id));
    for number=ids(:)'
        verifyGreaterThanOrEqual(testCase,number,1);
        verifyLessThanOrEqual(testCase,number,numel(registry.references));
        ref=registry.references(number);
        dois=as_cell(ref.dois);
        if isempty(dois)
            urls={ref.permalink};
        else
            urls=cellfun(@(doi)['https://doi.org/',doi],dois, ...
                'UniformOutput',false);
        end
        for j=1:numel(urls)
            verifyTrue(testCase,contains(source,lower(urls{j})), ...
                sprintf('Missing bibliography reference %d in %s',number,e.id));
        end
    end
end
end

function values=as_cell(value)
if isempty(value)
    values={};
elseif iscell(value)
    values=value(:).';
elseif isstring(value)
    values=cellstr(value);
elseif ischar(value)
    values={value};
else
    error('test_reference_registry:Type','Unexpected DOI list representation.');
end
end
