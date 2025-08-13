
% GetListOfFoldersInMasterFolder
% 
% When the script is run, the user is prompted to select a folder. The names of all the folders contained
% within that folder are written to an xlsx file named 'ListOfDatasetsInFolder' in the curent working directory.
% 
% Author: K Daniels
% Date: 19/04/19


ParentFolder = uigetdir('','Select Main Folder'); % Select the folder containing patients
ParentContent = dirPlus(ParentFolder, 'Depth', 0, 'ReturnDirs', true, 'PrependPath', false);
xlswrite ('ListOfDatasetsInFolder.xlsx', ParentContent);