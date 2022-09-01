function [varargout] = MOL_getTempPerSes(sesid,varargin)

for iArgin = 1:length(varargin)
    typeData = varargin{iArgin};
    
    %% Get the relevant data for each session individually:
    idx = false(size(typeData.session_ID));
%     for iSes = 1:length(sesid)
%         idx = idx | typeData.session_ID == sesid(iSes);
%     end
    
    for iSes = 1:length(sesid)
        idx = idx | strcmp(typeData.session_ID,sesid(iSes));
    end
    
    datafields = fieldnames(typeData);
    temptypeData = struct();
    for field = 1:length(datafields)
        temptypeData.(datafields{field}) = typeData.(datafields{field})(idx,:);
    end
    
    varargout{iArgin} = temptypeData;
end

end
