function [subfolders] = MOL_getsubfolders(folder)% get the folder contents
d = dir(folder);
% remove all files (isdir property is 0)
subfolders = d([d(:).isdir]==1);
% remove '.' and '..'
subfolders = subfolders(~ismember({subfolders(:).name},{'.','..'}));
subfolders = {subfolders.name};
% subfolders = fullfile(folder,subfolders);
end