% Moves or copies all files with a particular extension from one folder and
% its subfolders to another named folder. Change folder structure format if
% needed.

% Author: K Daniels
% Date: 11/06/1

% Choose top level folder to search for files within
start_path = 'C:\ViconDatabase';
topLevelFolder = uigetdir(start_path); % Generates GUI to allow top level directory to be selected
  %topLevelFolder = 'C:\ViconDataBase\Concussion\Baselinemissingtrials'; 

LIST = dirPlus(topLevelFolder, 'FileFilter', '\.csv'); 
% Returns cell array of all filepaths within folders and subfolders of topLevelFolder with FileFilter '.csv'

for n = 1:size(LIST,1)
    namepath = LIST{n} ;
    
    
  % Rename output file to include full path
    seps = find((namepath == '\')); % finds all backslashes in full file path
    path2 = namepath(seps(3)+1:end); % identifies the part of the file path that will be used in the new file name
    newfilename = strrep(path2,'\','_'); % Makes the new file name by changing backslashes to underscores
    newpath = fullfile('C:\ViconDataBase\Concussion\CSVs', newfilename) ; % creates the full new file path by appending the new file name to a destination folder name
  
  copyfile (LIST{n}, newpath) % copies all files in LIST to new folder defined in previous line, renaming each file to newfilename
end