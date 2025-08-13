function LabWeeklyNumbers_v2

% Function to extract Biomechanics Lab revelant weekly numbers from iMed-generated clinical report.

% Author: K Daniels, 
% Date: 01/03/17
% Updated to v2 by K Daniels 24/08/17 to include 'ISO Other' bookings in
% the ISO numbers.

% Function calls readindata (included in m file) authored by C Richter

%-----------------------------------------

% PROTOCOL FOR GENERATING WEEKLY NUMBERS SPREADSHEET

% 1) In iMed, select 'REPORT' tab' 
% 2) Select 'Clinical' then 'Clinic Appointments' from list on left of screen
% 3) Select 'Filter' (top right of page underneath the LOGOUT option)
% 4) Select From and To dates (Monday and Friday/Saturday of the relevant week)
% 5) Select 'Biomechanics' as the User Group. Leave Appointment Description and Status boxes as --Select--. 
% 6) Select '3D lab' in the Resource box
% 6) Click 'Export' and save where-ever you like. 
% 7) Run this code. Navigate to the exported report you have just saved when prompted.
% 8) Code will write an Excel file called 'labnumbers' to the folder in which you saved the iMed clinical report

%-----------------------------------------

% read in data
[filename,pathname] = uigetfile('*.xls','Select the Clinical Report file') ;% Navigate to and select clinical report, e.g. 'clinicdoc_report_9945218e4e7bdb306b1c92de7a8b24be.xls'
cd(pathname);

    data = readindata(filename);
    
 % determine UserGroups
    Type = '';
    for n = 1:size(data,1)
       var = data.Type{n}(5:end);
       var(strfind(var,' ')) = '';
       if ~any(strcmp(var,Type)) && ...
               ~sum(isnan(data.Type{n}))
           if ~isempty(strfind(var,'Groin3D')) || ...
                ~isempty(strfind(var,'ACL3D'))
            Type{size(Type,1)+1,1} = var;
           end
       end
       data.Type{n} = var;
    end
    
    
% extract relevant data and collate

    dateinfo = [];
    numinfo = [];
    
  idxday = find(data.SNO==1);
  for j = 1: size(idxday) 
      i = idxday(j);
   date = data.Date(i);
   
   if j<length(idxday)
   daydata = data(i:(idxday(j+1))-1,:) ; % get section of data for the day
   elseif i==idxday(end) % on last day of the selection
   daydata = data(i:end,:) ;    
   end
   
   numACLTest1_completed = sum((strcmp(daydata.Status, 'Completed') | strcmp(daydata.Status, 'Arrived'))& strcmp(daydata.Type, 'ACL3DTest1')); % number of ACLTest1 completed
   numACLTest2_completed = sum((strcmp(daydata.Status, 'Completed') | strcmp(daydata.Status, 'Arrived'))& strcmp(daydata.Type, 'ACL3DTest2')); % number of ACLTest2 completed
   numSport3D_completed = sum((strcmp(daydata.Status, 'Completed') | strcmp(daydata.Status, 'Arrived'))& strcmp(daydata.Type, 'Sport3D')); % number of Sport3D completed
   numACLTest1External_completed = sum((strcmp(daydata.Status, 'Completed') | strcmp(daydata.Status, 'Arrived'))& strcmp(daydata.Type, 'ACL3DTest1External')); % number of ACLTest1External completed
   numISO_completed = sum((strcmp(daydata.Status, 'Completed') | strcmp(daydata.Status, 'Arrived'))& (strcmp(daydata.Type, 'ISOTest-SurgeonReferral') | strcmp(daydata.Type, 'ISOother'))); % number of additional isos completed  
   numInitialGroin3D_completed = sum((strcmp(daydata.Status, 'Completed') | strcmp(daydata.Status, 'Arrived'))& strcmp(daydata.Type, 'InitialGroin3D')); % number of InitialGroin3D completed
   numGroin3Dretest_completed = sum((strcmp(daydata.Status, 'Completed') | strcmp(daydata.Status, 'Arrived'))& strcmp(daydata.Type, 'Groin3Dretest')); % number of Groin3Dretest completed
      
   numACLTest1_cancelled = sum((strcmp(daydata.Status, 'Cancelled') )& strcmp(daydata.Type, 'ACL3DTest1')); % number of ACLTest1 cancelled
   numACLTest2_cancelled = sum((strcmp(daydata.Status, 'Cancelled') )& strcmp(daydata.Type, 'ACL3DTest2')); % number of ACLTest2 cancelled
   numSport3D_cancelled = sum((strcmp(daydata.Status, 'Cancelled') )& strcmp(daydata.Type, 'Sport3D')); % number of Sport3D cancelled
   numACLTest1External_cancelled = sum((strcmp(daydata.Status, 'Cancelled') )& strcmp(daydata.Type, 'ACL3DTest1External')); % number of ACLTest1External cancelled
   numISO_cancelled = sum((strcmp(daydata.Status, 'Cancelled') )& (strcmp(daydata.Type, 'ISOTest-SurgeonReferral') | strcmp(daydata.Type, 'ISOother'))); % number of additional isos cancelled 
   numInitialGroin3D_cancelled = sum((strcmp(daydata.Status, 'Cancelled') )& strcmp(daydata.Type, 'InitialGroin3D')); % number of InitialGroin3D cancelled
   numGroin3Dretest_cancelled = sum((strcmp(daydata.Status, 'Cancelled') )& strcmp(daydata.Type, 'Groin3Dretest')); % number of Groin3Dretest cancelled
   
   numACLTest1_DNA = sum((strcmp(daydata.Status, 'DNA') )& strcmp(daydata.Type, 'ACL3DTest1')); % number of ACLTest1 DNA
   numACLTest2_DNA = sum((strcmp(daydata.Status, 'DNA') )& strcmp(daydata.Type, 'ACL3DTest2')); % number of ACLTest2 DNA
   numSport3D_DNA = sum((strcmp(daydata.Status, 'DNA') )& strcmp(daydata.Type, 'Sport3D')); % number of Sport3D DNA
   numACLTest1External_DNA = sum((strcmp(daydata.Status, 'DNA') )& strcmp(daydata.Type, 'ACL3DTest1External')); % number of ACLTest1External DNA
   numISO_DNA = sum((strcmp(daydata.Status, 'DNA') )& (strcmp(daydata.Type, 'ISOTest-SurgeonReferral') | strcmp(daydata.Type, 'ISOother'))); % number of additional isos DNA
   numInitialGroin3D_DNA = sum((strcmp(daydata.Status, 'DNA') )& strcmp(daydata.Type, 'InitialGroin3D')); % number of InitialGroin3D DNA
   numGroin3Dretest_DNA = sum((strcmp(daydata.Status, 'DNA') )& strcmp(daydata.Type, 'Groin3Dretest')); % number of Groin3Dretest DNA
   
   dateinfo = [dateinfo; date];
   numinfo = [numinfo; [numACLTest1_completed numACLTest2_completed numSport3D_completed numACLTest1External_completed...
       numISO_completed  numInitialGroin3D_completed numGroin3Dretest_completed...
       nan numACLTest1_cancelled numACLTest2_cancelled numSport3D_cancelled numACLTest1External_cancelled...
       numISO_cancelled  numInitialGroin3D_cancelled numGroin3Dretest_cancelled...
       nan numACLTest1_DNA numACLTest2_DNA numSport3D_DNA numACLTest1External_DNA...
       numISO_DNA  numInitialGroin3D_DNA numGroin3Dretest_DNA]   ] ;
   
  end
  
   dateinfo = dateinfo';
   numinfo = numinfo';
   total = sum(numinfo,2);
 totalcomp = sum(total(1:7));
 totalcanc = sum(total(9:15));
 totaldna = sum(total(17:23));
 subtotals = [totalcomp;totalcanc;totaldna];
 summarytitles = {'Total completed'; 'Total cancelled'; 'Total DNA'};
   
   headermain = {'Biomechanics Lab Weekly Numbers'};
   headerrow = dateinfo;
   headercol1 =  {'Completed'; nan;nan;nan;nan;nan;nan;nan; 'Cancelled'; nan;nan;nan;nan;nan;nan;nan; 'DNA'; nan;nan;nan;nan;nan;nan;nan   };
   headercol2 = {'ACLTest1'; 'ACLTest2'; 'Sport3D'; 'ACLTest1External'; 'ISO'; 'InitialGroin3D'; 'Groin3Dretest'; nan;... 
       'ACLTest1'; 'ACLTest2'; 'Sport3D'; 'ACLTest1External'; 'ISO'; 'InitialGroin3D'; 'Groin3Dretest'; nan;... 
       'ACLTest1'; 'ACLTest2'; 'Sport3D'; 'ACLTest1External'; 'ISO'; 'InitialGroin3D'; 'Groin3Dretest'};
    
  xlswrite('labnumbers.xlsx', headermain, 1, 'A1')
  xlswrite('labnumbers.xlsx', headerrow, 1, 'C3')
  xlswrite('labnumbers.xlsx', headercol1, 1, 'A4') 
  xlswrite('labnumbers.xlsx', headercol2, 1, 'B4') 
  xlswrite('labnumbers.xlsx', numinfo, 1, 'C4') 
    
  xlswrite('labnumbers.xlsx', {'TOTAL'}, 1, 'J3') 
  xlswrite('labnumbers.xlsx', total, 1, 'J4')    
    
 xlswrite('labnumbers.xlsx', summarytitles, 1, 'A29')      
 xlswrite('labnumbers.xlsx', subtotals, 1, 'B29')  
      
end




function data = readindata(file)

    [~,~,raw] = xlsread(file);
    data      = dataset;
    lables    = raw(11,:); 
    idx       = strcmp(raw(:,10),'Biomechanics');
    raw = raw(idx,:);

    % -------------------
    for n = 1:size(raw,2)
        
        isNumber = 1;
        for i = 2:size(raw,1)
            if ~isnumeric(raw{i,n})
                isNumber = 0;
                break
            end
        end

        VarName = lables{1,n};
        VarName(strfind(VarName,' ')) ='';
        VarName(strfind(VarName,'.')) ='';
        
        if strcmp(VarName,'AppointmentDescription')
            VarName = 'Type';
        end

        if isNumber
            data.(VarName) = cell2mat(raw(1:end,n));  
        else
            data.(VarName) = raw(1:end,n);
        end
        
    end
    
end         