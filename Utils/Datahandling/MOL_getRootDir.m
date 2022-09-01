function RootDir = MOL_getRootDir
%% Function to automatically locate the drive based on user pc:
hostname = strtrim(evalc('system(''hostname'');'));
if strcmp(hostname,'PCMatthijs')
    RootDir = 'E:';
elseif strcmp(hostname,'User-PC')
    RootDir = 'D:';
elseif any(strfind(hostname,'csn'))
    RootDir = '/data/moudelo1/';
else
    fprintf('Couldnt locate correct root directory')
end

end