
% DeleteC3DsInFolder

% Function deletes all C3D files in selected folder
% 
% Use: Use to delete processed data if second subject model accidentally
% attached or subject needs to be renamed.
%
% Author: K Daniels 
% Date: 02/03/17


    % --------------------------
    
%    % Detect all C3D files    
    
    origin = uigetdir('','Select the Session folder'); % Select the Session folder containing the c3d files you want to delete
    files = [origin filesep '*.c3d']; % Define file type to delete by extension name
    delete(files) % Delete files