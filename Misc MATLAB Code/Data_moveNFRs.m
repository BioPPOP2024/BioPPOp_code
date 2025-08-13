% Data_moveNFRs

% Script moves all patient folders with the string 'NFR' or 'Exclude'
% anywhere in the patient ENF file to alternative destination folders as
% defined below.
% Skips over any folder with multiple enf files at patient level and
% displays on screen.

% Author: K. Daniels
% Date: 25/03/19


originfolder = 'E:\02 LOWER LIMB PROCESSED DATA REPOSITORY\CLEANED collated'; % Origin folder
NFRdestinationfolder = 'E:\02 LOWER LIMB PROCESSED DATA REPOSITORY\NFR' ;  % Folder to move patient folder into if NFR
EXCLUDEdestinationfolder = 'E:\02 LOWER LIMB PROCESSED DATA REPOSITORY\Exclude' ;  % Folder to move patient folder into if excluded

folders = dirPlus(originfolder, 'Depth', 0, 'ReturnDirs', true); % list of patient folders in origin

for i = 1:length(folders) % loop through folders
    folders{i}
    
    % Check patient ENF
    patientEnf = dirPlus(folders{i}, 'Depth', 0, 'ReturnDirs', false, 'FileFilter', '.enf');
    
    if sum(size(patientEnf)) ~= 2 % Check that there is only one patient ENF
        disp('Wrong number of patient ENF files - check')
    continue
    end
    
    % Open patient ENF: 
              % open enffile
            fileID = fopen(char(patientEnf),'r');
          % read file line by line
            criteria = true; rep = 1; newENF = '';
            while criteria
                newENF{rep,1} = fgets(fileID);
                if newENF{rep} == -1
                    criteria = false;
                    newENF = newENF(1:rep-1,1);
                end
                if rep == 1
                    seperator = newENF{rep,1}(end-1:end);
                end
                    rep = rep + 1;
            end
          % close enffile                
            fclose(fileID);
            
 % Search whole ENF file for 'NFR' or 'exclude' (not case sensitive):
  if sum(contains(newENF, 'nfr', 'IgnoreCase',true)) > 0
      nfr = 1;
  else 
      nfr = 0;      
        if sum(contains(newENF, 'exclude', 'IgnoreCase',true)) > 0
        exclude = 1;
        else 
        exclude = 0;
        end      
  end
   
   % Move folder to NFR or Exclude folders if fulfil criteria:
  if nfr == 1
      status = movefile(folders{i}, NFRdestinationfolder);
  elseif exclude == 1
      status = movefile(folders{i}, EXCLUDEdestinationfolder);
  end
      
   
      
      
      
    
end

  
  
  