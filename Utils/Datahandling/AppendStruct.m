% AppendStruct appends trialdata, spike data or sessiondata to
% existing structs.
%    SOUT = AppendStruct(S, Sapp) appends in main struct S
%    all the structure fields, and their values,  that exist in
%    Sapp. If they that do not exist in S, they will be created and
%    aligned. If fields do not exist in Sapp they will be extended to align
%    with Sapp.

function sout = AppendStruct(s,sapp)
%% Do some checks on the input struct:

%% Make sure that fields to append are cell strings and not char arrays
sappfields = fieldnames(sapp);
for field = 1:length(sappfields)
    if ischar(sapp.(sappfields{field}))
        sapp.(sappfields{field}) = cellstr(sapp.(sappfields{field})); %Convert to cell if it was string
    end
end

%% Make sure that sessiondata fields to append do not contain multiple numeric elements
if isfield(s,'Mouse') || isfield(s,'mousename') || isfield(sapp,'Mouse') || isfield(sapp,'mousename') 
    for field = 1:length(sappfields)
        if isnumeric(sapp.(sappfields{field})) && ~all(size(sapp.mousename)==size(sapp.(sappfields{field}))) %numel(sapp.(sappfields{field}))> 1
            temp = sapp.(sappfields{field}); %Convert to cell if it was string
            sapp.(sappfields{field}) = cell(1,1); %Convert to cell if it was string
            sapp.(sappfields{field}){1} = temp;
        end
    end
end

%% Make sure that fields to append are not logical --> convert to numeric
sappfields = fieldnames(sapp);
for field = 1:length(sappfields)
    if islogical(sapp.(sappfields{field}))
        sapp.(sappfields{field}) = double(sapp.(sappfields{field}));
    end
end

%% Make sure that fields to append are not structs themselves
sappfields = fieldnames(sapp);
for field = 1:length(sappfields)
    if isstruct(sapp.(sappfields{field}))
        warning('Removing complete structure field')
        sapp = rmfield(sapp,sappfields{field});
    end
end

% %% Some parameters
% originalsize_s = unique(structfun(@length,s));
% originalsize_sapp = unique(structfun(@length,sapp));

%% Some parameters
originalsize_s = unique(structfun(@(x) size(x,1),s));
originalsize_sapp = unique(structfun(@(x) size(x,1),sapp));

%% Check that all fields of main struct s have the same number or rows
% if numel(unique(structfun(@length,s))) > 1
%     error('StructSize')
% end
if numel(unique(structfun(@(x) size(x,1),s))) > 1
    error('StructSize')
end

%% Check that all fields of struct to append s_app have the same number or rows
% if numel(unique(structfun(@length,sapp))) > 1
%     error('StructSize')
% end
if numel(unique(structfun(@(x) size(x,1),sapp))) > 1
    error('StructSize')
end

%% Make sure that fields to append are present in main struct
sappfields = fieldnames(sapp);
for field = 1:length(sappfields)
    if ~isfield(s,sappfields{field})
        if isnumeric(sapp.(sappfields{field}))
            s.(sappfields{field}) = NaN(originalsize_s,1);
        elseif islogical(sapp.(sappfields{field}))
            s.(sappfields{field}) = false(originalsize_s,1);
        elseif iscell(sapp.(sappfields{field}))
            s.(sappfields{field}) = cell(originalsize_s,1);
        elseif ischar(sapp.(sappfields{field}))
            sapp.(sappfields{field}) = cellstr(sapp.(sappfields{field})); %Convert to cell if it was string
            s.(sappfields{field})    = cell(originalsize_s,1);
        else     error('FieldType')
        end
    end
end

%% Make sure that all fields are of the same type:
for field = 1:length(sappfields)
    if isnumeric(sapp.(sappfields{field})) && iscell(s.(sappfields{field}))
        sapp.(sappfields{field}) = num2cell(sapp.(sappfields{field}));
    elseif iscell(sapp.(sappfields{field})) &&  isnumeric(s.(sappfields{field}))
        s.(sappfields{field}) = num2cell(s.(sappfields{field}));
    end
end
 
%% Append all fields
for field = 1:length(sappfields)
    s.(sappfields{field}) = [s.(sappfields{field}); sapp.(sappfields{field})];
end

%% Make sure that all fields of the main struct are aligned
sfields = fieldnames(s);
for field = 1:length(sfields)
    if ~isfield(sapp,sfields{field})
        if isnumeric(s.(sfields{field}))
            s.(sfields{field}) = [s.(sfields{field}); NaN(originalsize_sapp,1);];
        elseif iscell(s.(sfields{field}))
            s.(sfields{field}) = [s.(sfields{field}); cell(originalsize_sapp,1)];
        elseif islogical(s.(sfields{field}))
            s.(sfields{field}) = [s.(sfields{field}); false(originalsize_sapp,1)];
        else     error('FieldType')
        end
    end
end

%% Make sure that cell fields are not NaNs, but empty
for field = 1:length(sfields)
    if iscell(s.(sfields{field}))
        for iCell = 1:length(s.(sfields{field}))
            if isnan(s.(sfields{field}){iCell})
                s.(sfields{field}){iCell} = [];
            end
        end
    end
end

%% Check that all fields of main struct s have the same number or rows
% if numel(unique(structfun(@length,s))) > 1
%     error('StructSize')
% end
if numel(unique(structfun(@(x) size(x,1),s))) > 1
    error('StructSize')
end
sout = s;

end

