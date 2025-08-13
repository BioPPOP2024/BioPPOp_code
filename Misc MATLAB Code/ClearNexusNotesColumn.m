
% ClearNexusNotesColumn

% Function clears Notes from all trials in folder defined in line 14 and
% al its subfolders.
% 
% Use: Run on all data for full processing before giving out. 
%
% Author: K Daniels 
% Date: 02/03/17


    % --------------------------
    enftool   = EnfToolbox;
   % Detect all enf files
    origin = uigetdir('','Select the Session folder');
%     files = searchFolder4Files('I:\2016 DATA\TO CLEAN\March 2017') ; % get list of all files in folder and subfolders
    files = searchFolder4Files(origin) ; % get list of all files in folder and subfolders
   
        enftool   = EnfToolbox;
        DeleteIDX = zeros(size(files));
        for nfile = 1:size(files,1) % loop through all files

            if ~isempty(strfind(files{nfile,1},'.Trial')) && ~isempty(strfind(files{nfile,1},'.enf'))
                 DeleteIDX(nfile,:) = 1   ;
            end

        end
        files = files(DeleteIDX == 1); % List of trial ENF files
    
        for nfile = 1:size(files,1) % loop through all trial enf files
            enftool.changeNote(files{nfile,1},'')
        end
 
        
        
       