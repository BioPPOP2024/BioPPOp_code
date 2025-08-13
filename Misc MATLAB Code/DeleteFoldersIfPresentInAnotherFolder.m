% Written by Kat, 02/12/19.
% Function deletes any folder present in another named folder. Use to exclude patients from datasets to be processed if they have already been processed or flagged as NFRes

MainFolder = 'E:\WorkingFolder\ToProcess_Mar20'; % Folder with all the data in.
FolderToCheck = 'E:\02 LOWER LIMB PROCESSED DATA REPOSITORY\CLEANED collated' ; % Folder containing datasets that you want removed from MainFolder if present here.

cd(MainFolder)

MainContent = dirPlus(MainFolder, 'Depth', 0, 'ReturnDirs', true, 'PrependPath', false);
CheckContent = dirPlus(FolderToCheck, 'Depth', 0, 'ReturnDirs', true, 'PrependPath', false);

matchesID = contains(MainContent, CheckContent);

MainFull = dirPlus(MainFolder, 'Depth', 0, 'ReturnDirs', true, 'PrependPath', true); 
matchesList = MainFull(matchesID);

for i = 1:length(matchesList)
    toRem = char(matchesList(i));
    [status,msg] = rmdir(toRem, 's')
end