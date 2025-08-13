classdef EnfToolbox
%     
%     import ViconNexus
% vicon = ViconNexus()
%     
    properties
        % none
    end
    
    methods
        
        % find enf files --------------------------------------------------
        function EnfFile = GetEnfFiles(~,FILE)
            
          % transform to origin
            name  = FILE(max(strfind(FILE,filesep))+1:end-4);
            file  = FILE(1:max(strfind(FILE,filesep))-1);
            Session = file(max(strfind(file,filesep))+1:end);
            file  = file(1:max(strfind(file,filesep))-1);
            files = searchFolder4Files(file,'');
            
          % find enf files
            EnfFile = struct;
            for n = 1:size(files,1)
                
               lastSep  = max(strfind(files{n},filesep))+1;
               firstDot = min(strfind(files{n}(lastSep:end),'.'))+lastSep-2;
               currentfile = files{n}(lastSep:firstDot);
               
               if strcmp(currentfile,name) && ...
                       strcmp(files{n}(end-3:end),'.enf')  && ...
                       ~isempty(strfind(files{n},Session)) 
                EnfFile.TrialEnf = files{n};
               elseif ~isempty(strfind(files{n},'.Session')) && ...
                       strcmp(files{n}(end-3:end),'.enf')  && ...
                       ~isempty(strfind(files{n},Session))  
                EnfFile.SessionEnf = files{n};
               elseif ~isempty(strfind(files{n},'.Patient')) && ...
                       strcmp(files{n}(end-3:end),'.enf') 
                EnfFile.PatientEnf = files{n};
               end                
            end
            
            if sum(strcmp(fieldnames(EnfFile),'TrialEnf')) == 0
               EnfFile.TrialEnf = makeTrialEnf(FILE);
            end
            if sum(strcmp(fieldnames(EnfFile),'PatientEnf')) == 0
               EnfFile.PatientEnf = makePatientEnf(FILE);
            end
            
        end
        
        % change note section ---------------------------------------------
        function changeNote(~,enffile,term)
            
          % open enffile
            fileID = fopen(enffile,'r');
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
          % check if enf has already a note
            if sum(cell2mat(strfind(newENF,'NOTES='))) == 1
                count = 1;
                while count <= size(newENF,1)
                    if ~isempty(strfind(newENF{count,1},'NOTES='))
                        newENF{count,1} = ['NOTES=',term,seperator];
                    end
                    count = count+1;
                end
            elseif sum(cell2mat(strfind(newENF,'NOTES='))) == 0  
                newENF = addNoteToENFfile(enffile,newENF,term,seperator);
            elseif sum(cell2mat(strfind(newENF,'NOTES='))) >= 2
                count = 1; idx = zeros(size(newENF));
                while count <= size(newENF,1)
                    if ~isempty(strfind(newENF{count,1},'NOTES='))
                        idx(count) = 1;
                    end
                    count = count+1;
                end
                newENF = newENF(idx==0);
                newENF = addNoteToENFfile(enffile,newENF,term,seperator);
            end
            
            % check if note is behind patient info
            if ~isempty(strfind(enffile,'.Patient'))
                
                idx = zeros(size(newENF));
                DeleteNote = 'YES';
                for count = 1:size(newENF,1)
                    if ~isempty(strfind(newENF{count,1},'NOTES='))
                        idx(count) = 1;
                    elseif ~isempty(strfind(newENF{count,1},'[PATIENT_INFO]')) && ...
                            sum(idx) == 0
                        DeleteNote = 'NO';
                    end
                end
                
                if strcmp(DeleteNote,'YES')
                    newENF = newENF(idx == 0);
                    newENF = addNoteToENFfile(enffile,newENF,term,seperator);
                end
                
            end
            
          % overwrite enffile with new info                
            txt = fopen(enffile,'w');   
            for n = 1:size(newENF,1)
                fprintf(txt,'%s',newENF{n,:});
            end
            fclose(txt);
        end
        function changeDescription(~,enffile,term)
          % open enffile
            fileID = fopen(enffile,'r');
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
          % check if enf has already a note
            if sum(cell2mat(strfind(newENF,'DESCRIPTION='))) == 1
                count = 1;
                while count <= size(newENF,1)
                    if ~isempty(strfind(newENF{count,1},'DESCRIPTION='))
                        newENF{count,1} = ['DESCRIPTION=',term,seperator];
                    end
                    count = count+1;
                end
            elseif sum(cell2mat(strfind(newENF,'DESCRIPTION='))) == 0  
                
                newENF = addDescribtionToENFfile(enffile,newENF,term,seperator);
                
            elseif sum(cell2mat(strfind(newENF,'DESCRIPTION='))) >= 2
                
                count = 1; idx = 0;
                while count <= size(newENF,1)
                    if ~isempty(strfind(newENF{count,1},'DESCRIPTION='))
                        if idx == 0
                            newENF = newENF(1:count-1,:);
                        end
                    end
                    count = count+1;
                end
                newENF = addDescribtionToENFfile(enffile,newENF,term,seperator);
                
            end
            
            % check if note is behind patient info
            if ~isempty(strfind(enffile,'.Patient'))
                
                idx = zeros(size(newENF));
                DeleteNote = 'YES';
                for count = 1:size(newENF,1)
                    if ~isempty(strfind(newENF{count,1},'DESCRIPTION='))
                        idx(count) = 1;
                    elseif ~isempty(strfind(newENF{count,1},'[PATIENT_INFO]')) && ...
                            sum(idx) == 0
                        DeleteNote = 'NO';
                    end
                end
                
                if strcmp(DeleteNote,'YES')
                    newENF = newENF(idx == 0);
                    newENF = addDescribtionToENFfile(enffile,newENF,term,seperator);
                end
                
            end
          % overwrite enffile with new info                
            txt = fopen(enffile,'w');   
            for n = 1:size(newENF,1)
                fprintf(txt,'%s',newENF{n,:});
            end
            fclose(txt);
        end
        
        % check note section ----------------------------------------------
        function [criteria,term] = GetNote(~,file,screenterm,field)

          % open file
            fileID = fopen(file,'r');
          % read file line by line
            income = textscan(fileID,'%s','delimiter',',');
            income = income{1,1};
          % detect note section and the describtion
            fclose(fileID);

          % look for note section
            term = '';
            for n = 1:size(income,1)
                if ~isempty(strfind(income{n,1},field))
                    term = income{n,1};
                    term = term(size(field,2):end);
                end
            end

          % see if criteria is met ...
            if iscell(screenterm)
                for n = 1:size(screenterm,2)
                  if ~isempty(strfind(term,screenterm{n}))
                      criteria = true;
                      break
                  else
                      criteria = false; 
                  end
                end
            else
              if ~isempty(strfind(term,screenterm))
                  criteria = true;
              else
                  criteria = false; 
              end 
            end
            if strcmp(screenterm,'none')
                criteria = true;
            end
            
        end
        
    end
    
end

% help functions ---------------------------------------------------
    % get all files of folder and subfolder ------------------------
function files = searchFolder4Files(origin,files)
    folder = dir(origin);
    for n = 3:size(folder,1)
            suborigin = [origin,filesep,folder(n).name];
        switch isdir(suborigin)
          case 1
            files = searchFolder4Files(suborigin,files);
          case 0
            if ~isempty(strfind(suborigin,'.enf'))
             files{size(files,1)+1,1} = suborigin;
            end
        end       
    end
end
    % add note if not existing -------------------------------------
function newENF = addNoteToENFfile(enffile,newENF,term,seperator)

    if ~isempty(strfind(enffile,'.Trial'))
    % patient and session enf need to be treaded differently
        for n=1:size(newENF,1) 
            if strfind(newENF{n},'CREATIONDATEANDTIME')
                part1 = newENF(1:n);
                part2 = newENF(n+1:end,1);
                newENF = [part1;'NOTES=',term,seperator;part2];
                break
            end
            if size(newENF,1) == n 
                part1 = newENF(1:n-2);
                part2 = newENF(n-1:end,1);
                newENF = [part1;'NOTES=',term,seperator;part2];
            end
        end    
    elseif ~isempty(strfind(enffile,'.Session'))

        for n=1:size(newENF,1) 
            if strfind(newENF{n},'CREATIONDATEANDTIME')
                newENF{n+1} = ['NOTES=',term,seperator];
                newENF(n+2:size(newENF,1)+1) = newENF(n+1:end);
                break
            end
        end

    elseif ~isempty(strfind(enffile,'.Patient'))
        
        for n=1:size(newENF,1) 
            if strfind(newENF{n},'[PATIENT_INFO]')
                part1 = newENF(1:n);
                part2 = newENF(n+1:end,1);
                newENF = [part1;'NOTES=',term,seperator;part2];
                break
            end
            if size(newENF,1) == n 
                newENF = [newENF;'[PATIENT_INFO]',seperator;'NOTES=',term,seperator];
            end
        end

    else
        error('not defined')
    end

end
function newENF = addDescribtionToENFfile(enffile,newENF,term,seperator)

    if ~isempty(strfind(enffile,'.Trial')) || ...
        ~isempty(strfind(enffile,'.Session')) || ...
        ~isempty(strfind(enffile,'.Patient'))
        
         for n=1:size(newENF,1) 
            if strcmp(newENF{n},['[PATIENT_INFO]',seperator])
                part1 = newENF(1:n);
                part2 = newENF(n+1:end,1);
                newENF = [part1;'DESCRIPTION=',term,seperator;part2];
                break
            end
            if size(newENF,1) == n 
                newENF = [newENF;...
                          '[PATIENT_INFO]',seperator;...
                          'DESCRIPTION=',term,seperator]; %#ok<AGROW>
            end
        end
        
    else
       
        error('not defined')
       
    end
end
    % make new enf file --------------------------------------------
function Enf = makeTrialEnf(file)

    origin = which('EnfToolbox');
    origin = origin(1:max(strfind(origin,filesep)));
    origin = [origin 'ExampleFilesEnfToolBox'];
    if exist([pwd filesep 'Example.Trial01.enf'],'file') == 2
        enfFile = [pwd filesep 'Example.Trial01.enf'];
    elseif exist([origin filesep 'Example.Trial01.enf'],'file') == 2
        enfFile = [origin filesep 'Example.Trial01.enf'];
    else
       warndlg('No Trial Enf File');
       return 
    end
    SName = file(max(strfind(file,filesep)+1):end);
    
  % open enffile
    fileID = fopen(enfFile,'r');
  % read file line by line
    criteria = true; rep = 1; newENF = '';
    while criteria
        newENF{rep,1} = fgets(fileID);
        if ~isempty(strfind(newENF{rep},'SUBJECTS=Example'))
            income = newENF{rep,1};
            newENF{rep,1} = ['SUBJECTS=',SName,income(end-1:end)];
        elseif newENF{rep} == -1
            criteria = false;
            newENF = newENF(1:rep-1,1);
        end

        rep = rep + 1;
    end
  % close enffile                
    fclose(fileID);
        
  % write enffile with new info   
    Enf = [file(1:end-3),'Trial01.enf'];
    txt = fopen(Enf,'w');   
    for n = 1:size(newENF,1)
        fprintf(txt,'%s',newENF{n,:});
    end
    fclose(txt);
    
end
function Enf = makePatientEnf(file)

    origin = which('EnfToolbox');
    origin = origin(1:max(strfind(origin,filesep)));
    origin = [origin 'ExampleFilesEnfToolBox'];
    if exist([pwd filesep 'Example.Patient 1.enf'],'file') == 2
        enfFile = [pwd filesep 'Example.Trial01.enf'];
    elseif exist([origin filesep 'Example.Patient 1.enf'],'file') == 2
        enfFile = [origin filesep 'Example.Patient 1.enf'];
    else
       warndlg('No Patient Enf File');
       return 
    end
    SName = file(max(strfind(file,filesep)+1):end);
    
  % open enffile
    fileID = fopen(enfFile,'r');
  % read file line by line
    criteria = true; rep = 1; newENF = '';
    while criteria
        newENF{rep,1} = fgets(fileID);
        if ~isempty(strfind(newENF{rep},'NAME='))
            income = newENF{rep,1};
            newENF{rep,1} = ['SUBJECTS=',SName,income(end-1:end)];
        elseif ~isempty(strfind(newENF{rep},'SESSION='))
            income = newENF{rep,1};
            newENF{rep,1} = ['SESSION=',SName,income(end-1:end)];
        elseif newENF{rep} == -1
            criteria = false;
            newENF = newENF(1:rep-1,1);
        end

        rep = rep + 1;
    end
  % close enffile                
    fclose(fileID);
        
  % write enffile with new info   
    Enf = [file(1:end-3),'Patient 1.enf'];
    txt = fopen(Enf,'w');   
    for n = 1:size(newENF,1)
        fprintf(txt,'%s',newENF{n,:});
    end
    fclose(txt);
    
end