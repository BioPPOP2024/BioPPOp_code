% Transfers and copies data of defined file type.

% Load spreadsheet containing info regarding target files
clear;

%file = 'ACL_Master 23 May_withSurveyMonkeyDataAndCorrectedQuestionnaires.xlsx';
file = 'F:\ProcessedStudyData\ACL Reruptures\FemaleControls_April20\Female Controls';
datadirectory = 'F:\01 RAW DATA BACKUP\ACLGROIN\AllData';
originfolder = datadirectory;
cd(datadirectory)
%originfolder = cd;
destinationfolder = 'F:\ProcessedStudyData\ACL Reruptures\FemaleControls_April20\FemaleControls_ToBeProcessed' ;  % folder to copy data into


    %[~,~,raw] = xlsread(file);  
    [~,~,raw] = xlsread(file, 'ToBeProcessed');  
    
    List = dataset;
    for n = 1:size(raw,2) 
        
        % check if number 
        IsNumber = 1;
        for i = 2:size(raw,1)
            if ~isnumeric(raw{i,n})
                IsNumber = 0;
            end
        end
        
        % export
        colheader = raw{1,n};
        colheader = strrep(colheader,' ','');
        if IsNumber
          %  List.(raw{1,n}) = cell2mat(raw(2:end,n));
            List.(colheader) = cell2mat(raw(2:end,n));
        else
          %  List.(raw{1,n}) = raw(2:end,n);
            List.(colheader) = raw(2:end,n);
        end
    end
    
    % Create additional columns to indicate whether data are present
    List2 = List;
    nrow = size(List,1);
    List2.Test1_Data = zeros(nrow, 1);
    List2.Test2_Data = zeros(nrow, 1);
    
    
    for i = 1:size(List2,1)
        try
       %ssc = List2.PatientID_PreOp(i);
       ssc = List2.PatientID(i);
       % Strip letters and leading zeros to get SSC as per Vicon data
       % ACL naming conventions and add leading whitespace:
       ssc = char(ssc); ssc = ssc(end-6:end); ssc(regexp(ssc(1:2),'0'))=[]; ssc = [' ', ssc];        
       % Get first initial as additional check:
       initial = List2.PatientName(i); initial = char(initial); initial = initial(:,1);
        
       % Search for patient datasets in current directory
       dinfo = dir();
       dinfo(ismember({dinfo.name}, {'.', '..'})) = [];  % remove '.' and '..'
       dinfo = struct2table(dinfo); % convert to table
       
       IndexC = strfind(dinfo.name, ssc)  ;      
       Index = find(not(cellfun('isempty', IndexC))); % returns index to row(s) of dinfo.name
       if isempty(Index); continue % skip to next dataset if no matching SSC found
       end
       
       % Check whether first initial matches:    
       checkInitial = strncmpi(dinfo.name(Index), initial, 1); % returns 1 for match, 0 for no match
       
       % Check whether revision or contra:
       checkRepeat = strcmpi(dinfo.name(Index), 'Revision') | strcmpi(dinfo.name(Index), 'Contra') ; 
       
      % Check whether retest: 
       checkRetest = strfind(dinfo.name(Index), 'etest') ; 
       checkRetest = not(cellfun('isempty', checkRetest));
        catch
            continue
        end
        

        c = 0;
       for k = Index' 
           c = c+1;
           %if checkInitial(c,:) == 1 && checkRepeat(c,:) == 0 && checkRetest(c,:) == 1; % only proceed if first initial matches and not a contra/revision
           %if  checkInitial(c,:) == 1 && checkRepeat(c,:) == 0 && checkRetest(c,:) == 0 % only proceed if first initial matches, not a contra/revision and not a retest
           if  checkInitial(c,:) == 1 && checkRepeat(c,:) == 0 && checkRetest(c,:) == 1 % only proceed if first initial matches, not a contra/revision and a retest
               dinfo.name(k)
               % Copy dataset over
               patientfolder = fullfile(originfolder, dinfo.name(k));
               dest = fullfile(destinationfolder, dinfo.name(k));
if 1==1              
               copyfile (patientfolder{1}, dest{1})
               
%                try
%                % Delete unwanted file types
%                enflist = dirPlus(dest{1}, 'FileFilter', '.enf');
%                delete(enflist{:,:}); 
%                systemlist = dirPlus(dest{1}, 'FileFilter', '.system');
%                delete(systemlist{:,:}); 
%                historylist = dirPlus(dest{1}, 'FileFilter', '.history');
%                delete(historylist{:,:}); 
%                mplist = dirPlus(dest{1}, 'FileFilter', '.mp');
%                delete(mplist{:,:});
%                vsklist = dirPlus(dest{1}, 'FileFilter', '.vsk');
%                delete(vsklist{:,:}); 
%                x1dlist = dirPlus(dest{1}, 'FileFilter', '.x1d');
%                delete(x1dlist{:,:}); 
%                x2dlist = dirPlus(dest{1}, 'FileFilter', '.x2d');
%                delete(x2dlist{:,:}); 
%                xcplist = dirPlus(dest{1}, 'FileFilter', '.xcp');
%                delete(xcplist{:,:}); 
%                catch
%                     warning('Not all expected file types present: some file types may not have been deleted')
%                end
%                
%                try
%                mkrlist = dirPlus(dest{1}, 'FileFilter', '.mkr');
%                delete(mkrlist{:,:}); 
%                catch
%                     warning('No mkr files present')
%                end
%                try
%                digdevlist = dirPlus(dest{1}, 'FileFilter', '.xml');
%                delete(digdevlist{:,:}); 
%                catch
%                     warning('No xml files present')
%                end
%                
               try
               avilist = dirPlus(dest{1}, 'FileFilter', '.avi');
               delete(avilist{:,:}); 
               catch
                    warning('No avi files present')
               end
               
               

 %              try
%                c3dlist = dirPlus(dest{1}, 'FileFilter', '.c3d'); % Gets list of all c3d files in folder
%                  %IndexEx1 = strfind(c3dlist, '\CMJ');
%                  IndexEx1 = strfind(c3dlist, 'Cut');
%                  IndexEx = find(not(cellfun('isempty', IndexEx1))); % Get indices to exercise                  
%                  IndexSt1 = strfind(c3dlist, '\Static');
%                  IndexSt = find(not(cellfun('isempty', IndexSt1))); % Get indices to Static Trial      
%                  
%                  IndexChosenTrials = [IndexEx; IndexSt];
%                  alltrials = (1:length(c3dlist))';                  
%                  IndexDeleteTrials = setdiff(alltrials, IndexChosenTrials); 
% 
%                unwantedc3dlist = c3dlist(IndexDeleteTrials, :); % List of c3d trials to delete
%                delete(unwantedc3dlist{:,:})
%                catch
%                    warning('Problem using function - some c3d files from unwanted exercises may not have been deleted')
%              
%                end

end               
               
               if checkRetest(c,:) == 0 % Write 1 into relevant column (Test1/Test2) to say that dataset was found, checked and copied over
                   List2.Test1_Data(i,:) = 1;
               elseif checkRetest(c,:) == 1
                   List2.Test2_Data(i,:) = 1;
               end
           else % If first initial does not match or dataset is contra/revision, write 9 into relevant List column
               if checkRetest(c,:) == 0
                   List2.Test1_Data(i,:) = 9;
               elseif checkRetest(c,:) == 1
                   List2.Test2_Data(i,:) = 9;
               end
               
           end
        end           
     
      
    end
    
     % Write new list to Excel spreadsheet
     warning('off','MATLAB:xlswrite:AddSheet')
  %    xlswrite('Hams091019.xlsx',List2,'Sheet1','A1');
    export(List2, 'file', 'FemaleControlsToBeProcessedList.txt')
    
    
    
    

