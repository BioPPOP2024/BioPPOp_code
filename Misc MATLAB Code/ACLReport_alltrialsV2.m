function ACLReport_alltrialsV2(origin) % ----------------------------------------------
% -------------------------------------------------------------------------

% This m-file lets the user select a folder. The m-file will search only 
% within the selected folder for c3d-files that has a specific term in  
% their note or description section. Every file that has been identified 
% and previously defined (see below) will then be used to generate an 
% ACLReport as requested by Marit Undheim

% -------------------------------------- Chris Richter, DCU, 16th Sep. 2014
% ---------------------------------------------- mr.chris.richter@gmail.com
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% changes 08/10/14 --- MTP has been added as requested by Marit Undheim ---
% -------------------- Format has changed slightly ------------------------
% -------------------------------------------------------------------------

% changes 17/11/14 --- fixed bug in Code: Crashed when not all exercises --
% -------------------- where found ----------------------------------------
% -------------------------------------------------------------------------

% changes 22/06/15 --- fixed bug in Code: Code can now handle subject -----
% -------------------- folder inputs --------------------------------------
% -------------------------------------------------------------------------

    if nargin<1
        origin = uigetdir;
    end
    
    close all

    Porigin = which('ACLReport_alltrialsV2');
    cd(Porigin(1:(max(strfind(Porigin,filesep)))-1))
    
    % please specifiy the term for file that have to be screened
        screenterm = 'U'; field = 'DESCRIPTION=';
        
% -------------------------------------------------------------------------       

    fprintf(' --- PROGRAM STARTING --- \n\n');    
    h = waitbar(0,'Please wait ... (searching for files)');
    
    % get all the c3d files within directory
        files = searchFolder4Files(origin); 
     
    waitbar(0,h,'Please wait ... (screening files)')
    
    % data processing -----------------------------------------------------
        
        measures     = dataset;
        nameIDX      = (strfind(origin,filesep));
        subject_name = origin(nameIDX(end-1)+1:end);
        subject_name(strfind(subject_name,filesep)) = '-';
        
    % find model file
    
    %find BW
     IDX = zeros(size(files,1),1);
    for n = 1:size(files)
        if ~isempty(strfind(files{n},'.mp'))
          modelfile = files{n};
          BW = getBodyWeight(modelfile);
          IDX(n,1) = 1;
        end
    end
        files = files(IDX==0);
        
    for sub = 1:size(files,1) 
        
        waitbar(sub / size(files,1))
        
        % check for c note in the files 
         %   useQuest = CheckNote(files{sub},...
%                                  screenterm,...
%                                  field);
        useQuest = 1; %% KD edit 08/01/16 to take all files
        if useQuest
            
            % this function extracts data all the c3d files detected...
                [data,info,type] = getDataANDstrucData(files{sub});

            if isempty(info)
                
                measures = data2discretepoints(data,type,measures,files{sub},BW);
                
                file = files{sub};
                file = file(size(origin,2)+2:end);
                file(strfind(file,filesep)) = '-';
                fprintf([subject_name,'-',file,'-',num2str(sub),'\n'])
                
            else
                
                if ~strcmp(info,'exercise no recognized')
                info(strfind(info,filesep)) = '-';
                fprintf(2,['\t',subject_name,info,'-',num2str(sub),'\n']);
                end
                
            end

        end
        
    end

    waitbar(1,h,'Please wait ... (writing report)')
    
    % bring into excel file -----------------------------------------------
    
    	makeExcelSheet(measures,origin)
    
    % ---------------------------------------------------------------------
    
    fprintf(' --- PROGRAM FINISHED --- \n');
    
    delete(h)
    
end

function files = searchFolder4Files(origin)

 folder = dir(origin);
 files  = '';
 
    for n = 3:size(folder,1)
        suborigin = [origin,filesep,folder(n).name];
        switch isdir(suborigin)
         case 1
             files = getintofolder(suborigin,files);
         case 0
             if ~isempty(strfind(suborigin,'.c3d')) || ...
                     ~isempty(strfind(suborigin,'.mp'))
             files{size(files,1)+1,1} = suborigin;
             end
        end       
    end
    
end

    function files = getintofolder(origin,files)

        folder = dir(origin);
        for n = 3:size(folder,1)
            suborigin = [origin,filesep,folder(n).name];
            switch isdir(suborigin)
             case 1
                 files = getintofolder(suborigin,files);
             case 0
                 if ~isempty(strfind(suborigin,'.c3d')) || ...
                     ~isempty(strfind(suborigin,'.mp'))
                 files{size(files,1)+1,1} = suborigin;    
                 end
            end       
        end

    end
    
function criteria = CheckNote(file,screenterm,field)

   % find enf.-
    origin   = dir(file(1:max(strfind(file,filesep))-1));
    filename = file(max(strfind(file,filesep))+1:end-4);
    
    for n = 1:size(origin,1)
       if ~isempty(strfind(origin(n).name,filename)) && ...
             ~isempty(strfind(origin(n).name,'.enf'))   

            count = n;
            
       end
    end

   % create enf file name
   enfFile = [file(1:max(strfind(file,filesep))),origin(count).name];
   % open file
   fileID = fopen(enfFile,'r');
   % read file line by line
   income = textscan(fileID,'%s','delimiter',',');
   income = income{1,1};
   % detect note section and the describtion
   fclose(fileID);
   
   term       = '';
   
   for n = 1:size(income,1)
      if ~isempty(strfind(income{n,1},field))
          term = income{n,1};
          term = term(size(field,2):end);
      end
   end
   
   % see if criteria is met ...
   if ~isempty(strfind(term,screenterm))
    criteria = true;
   else
    criteria = false; 
   end
   
end

function [data,info,type] = getDataANDstrucData(file)

    info = [];
    idx  = strfind(file,filesep);

    % set up variable for further processing
    if  ~isempty(strfind(file(idx(end)+1:end),'SL')) && ...
            ~isempty(strfind(file(idx(end)+1:end),'CMJ'))
        type = 'SL_CMJ';
    elseif ~isempty(strfind(file(idx(end)+1:end),'CMJ'))
        type = 'DL_CMJ';
    elseif~isempty(strfind(file(idx(end)+1:end),'SL')) && ...
            ~isempty(strfind(file(idx(end)+1:end),'Hop'))
        type = 'SL_Hop';
    elseif ~isempty(strfind(file(idx(end)+1:end),'DL')) && ...
            ~isempty(strfind(file(idx(end)+1:end),'DJ'))
        type = 'DL_DJ';
    elseif ~isempty(strfind(file(idx(end)+1:end),'SL')) && ...
            ~isempty(strfind(file(idx(end)+1:end),'DJ'))
        type = 'SL_DJ';
    elseif ~isempty(strfind(file(idx(end)+1:end),'SL')) && ...
            ~isempty(strfind(file(idx(end)+1:end),'DJ'))
        type = 'SL_DJ';
    elseif ~isempty(strfind(file(idx(end)+1:end),'Bilat MTP'))
        type = 'DL_MTP';
    elseif ~isempty(strfind(file(idx(end)+1:end),'Uni MTP'))
        type = 'SL_MTP';
    else
        type = 'not defined'; data = [];info = 'exercise no recognized';
    end
    
    if isempty(info)
        % define variables to look at
            if ~isempty(strfind(file(idx(end)+1:end),'Left'))
                side = 'L';
            elseif ~isempty(strfind(file(idx(end)+1:end),'Right'))
                side = 'R';
            else
                side = 'B';
            end

        if strcmp(type,'SL_CMJ') || ...
                strcmp(type,'DL_CMJ') || ...
                strcmp(type,'DL_DJ') || ...
                strcmp(type,'SL_DJ') || ...
                strcmp(type,'DL_MTP') || ...
                strcmp(type,'SL_MTP')

                 VarOfInt  = {'Force_Fz1','Force_Fz2'};

        elseif strcmp(type,'SL_Hop')

                 VarOfInt  = {'TOE','Force_Fz1','Force_Fz2'};

        else
                 error('Exercise is has no variables of interest defined')
        end

        % load data within incoming file

           c3dfile = btkReadAcquisition(file);
           marker  = btkGetMarkers(c3dfile);

           % get raw forces....
            analogs = btkGetAnalogs(c3dfile);
            ratio = btkGetAnalogSampleNumberPerFrame(c3dfile);
            analogsDownsampled = [];
            labels = fieldnames(analogs);
            for i = 1:btkGetAnalogNumber(c3dfile)
            analogsDownsampled.(labels{i}) = analogs.(labels{i})(1:ratio:end);
            end

       DATA = mergestruct(marker,analogsDownsampled);
       
       btkCloseAcquisition(c3dfile)

       % check for variables of interest in marker

        data   = dataset;
        info   = '';
        labels = fieldnames(DATA);

        for n_mar = 1:size(VarOfInt,2)

            if strcmp(VarOfInt{n_mar},'TOE')
                VARIABLE = [side,VarOfInt{n_mar}];
            else
                VARIABLE = VarOfInt{n_mar};
            end

            if sum(cell2mat(strfind(labels,VARIABLE))) >= 1
                data.(VarOfInt{n_mar}) = DATA.(VARIABLE);
            else
                if strcmp(type,'Static')
                data.(VarOfInt{n_mar}) = DATA.(labels{1}) + NaN;   
                else
                info = [info,'-',file(idx(end-1)+1:end),'-',VARIABLE,' missing, ']; %#ok<AGROW>
                end
            end
        end
    end
end

function measures = data2discretepoints(data,type,measures,file,BW)

    name = file(max(strfind(file,filesep))+1:end-4);

    if ~isempty(strfind(name,'Left'))
        side = 'L';
    elseif ~isempty(strfind(name,'Right'))
        side = 'R';
    else
        side = 'B';
    end
    
    % -----

    sample = dataset;

    % ---------------------------------------------------------------------

    if strcmp(type,'DL_CMJ') || strcmp(type,'SL_CMJ')
        
       if strcmp(side,'L')
          GRF = -data.Force_Fz2;  
       elseif strcmp(side,'R')
          GRF = -data.Force_Fz1;  
       else
          GRF = -(data.Force_Fz1 + data.Force_Fz2); 
       end
       
        takeoff   = find(GRF < 5,1,'first');%+1; edit by K Daniels, 17/06/15
    [~,threshold] = max(GRF(takeoff:end));
        threshold = threshold + takeoff - 1;
        impact    = find(GRF(1:threshold) < 5,1,'last');

        duration = (impact-takeoff)*(1/200);

        velo = (9.81*duration) / 2;
        jump_height = (velo^2 / (2*9.81))*100;
        
        % add variables to dataset 
        DP = jump_height;
        
        ds = dataset;
        ds.trial{1}    = name;        
        ds.task{1}     = type;
        ds.side{1}     = side;          
        ds.measures{1} = {'jump height'};
        ds.VarIDX      = 1;
        
        sample = [sample ds];

            sample = [sample dataset(DP)];
            sample.Properties.VarNames{6} = 'value';
            
        measures = [measures;sample];
            
    elseif strcmp(type,'DL_DJ') || strcmp(type,'SL_DJ')
        
       if strcmp(side,'L')
          GRF = -data.Force_Fz2;  
       elseif strcmp(side,'R')
          GRF = -data.Force_Fz1; 
       else
          GRF = -(data.Force_Fz1 + data.Force_Fz2); 
       end
       
        threshold = find(GRF > max(GRF)*.4,1,'first');        
        start     = find(GRF(1:threshold) < 5,1,'last');
        takeoff   = find(GRF(threshold:end) < 5,1,'first')+threshold-1;
    [~,threshold] = max(GRF(takeoff:end));
        threshold = threshold + takeoff - 1;
        impact    = find(GRF(1:threshold) < 5,1,'last')+1;

        duration = (impact-takeoff)*(1/200);

        velo = (9.81*duration) / 2;
        
        if isempty(impact) || isempty(takeoff)
           error('Pleasse check Vicon/Code (threshold) file as no start or takeoff was identified') 
        end
        
        jump_height = (velo^2 / (2*9.81))*100;
        contacttime = (takeoff-start)*(1/200);
        rsi = (jump_height./100)/contacttime ;
      
         
        % add variables to dataset 
        DP = [jump_height;contacttime;rsi];
            
        ds = dataset;
        ds.trial    = {name;name;name};        
        ds.task     = {type;type;type};
        ds.side     = {side;side;side};          
        ds.measures = {'jump height';'contact time';'RSI'};
        ds.VarIDX   = [1;2;3];
        
        sample = [sample ds];
            
        for n = 1:size(DP,2)

            sample = [sample dataset(DP(:,n))]; %#ok<AGROW>
            sample.Properties.VarNames{n+5} = 'value';

        end
            
        measures = [measures;sample];
            
    elseif strcmp(type,'SL_Hop')
        
        GRF = -(data.Force_Fz1 + data.Force_Fz2); 
        
        [~,peak] = max(GRF);
        
        toe = data.TOE(:,2) - min(data.TOE(:,2));
    
        distance = (mean(toe(1:peak))-mean(toe(end-20:end)))/10;
         
        % add variables to dataset 
        DP = distance;
            
        ds = dataset;
        ds.trial{1}    = name;        
        ds.task{1}     = type;
        ds.side{1}     = side;          
        ds.measures{1} = {'distance'};
        ds.VarIDX      = 1;
        
        sample = [sample ds];

            sample = [sample dataset(DP)];
            sample.Properties.VarNames{6} = 'value';
            
        measures = [measures;sample];
        
    elseif strcmp(type,'DL_MTP') || strcmp(type,'SL_MTP')
    
        GRF = -(data.Force_Fz1 + data.Force_Fz2); 
        
        if isempty(BW)
            peak = NaN;
        else
            peak = max(GRF) / BW;
        end
        
        % add variables to dataset 
        DP = peak;
            
        ds = dataset;
        ds.trial{1}    = name;        
        ds.task{1}     = type;
        ds.side{1}     = side;          
        ds.measures{1} = {'Peak Force'};
        ds.VarIDX      = 1;
        
        sample = [sample ds];

            sample = [sample dataset(DP)];
            sample.Properties.VarNames{6} = 'value';
            
        measures = [measures;sample];
        
    end       
end

function makeExcelSheet(measures,origin)

    warning('off','MATLAB:xlswrite:AddSheet')

    %make an IDX
    IDX = zeros(size(measures,1),2);

    for n = 1:size(measures,1)
       switch measures.task{n}
           case 'SL_CMJ';  IDX(n,1) = 1;
           case 'DL_CMJ';  IDX(n,1) = 2;
           case 'SL_Hop';  IDX(n,1) = 3;
           case 'DL_DJ';   IDX(n,1) = 4;
           case 'SL_DJ';   IDX(n,1) = 5;
           case 'DL_MTP';  IDX(n,1) = 6;
           case 'SL_MTP';  IDX(n,1) = 7;
           otherwise
               error(['Exercise not defined: ',measures.task(n)])
       end
       switch measures.side{n}
           case 'L';      IDX(n,2) = 1;
           case 'R';      IDX(n,2) = 2;
           case 'B';      IDX(n,2) = 3;
       end
    end
    
    % bring the variables in output format
    
        nameIDX = (strfind(origin,filesep));
        subject_name = origin(nameIDX(end-1)+1:end);
        subject_name(strfind(subject_name,filesep)) = '-';
       
        output = '';
        
    for nex = 1:7
        
        switch nex
            case 1; type = 'SL_CMJ';
            case 2; type = 'DL_CMJ';
            case 3; type = 'SL_Hop';
            case 4; type = 'DL_DJ';
            case 5; type = 'SL_DJ';
            case 6; type = 'DL_MTP';
            case 7; type = 'SL_MTP';
        end
        
        if isempty(strfind(type,'DL_'))
            
            Loutput = '';
        % left ...
            index = ((IDX(:,1)==nex) + (IDX(:,2)==1)) == 2;
            Loutput = data2string(Loutput,measures(index,:),['L ',type]);

            Routput = '';
        % right ...
            index = ((IDX(:,1)==nex) + (IDX(:,2)==2)) == 2;
            Routput = data2string(Routput,measures(index,:),['R ',type]);
            
            if size(Loutput,1) > size(Routput,1)
              filler = size(Loutput,1) - size(Routput,1);
              Routput = [Routput;Routput(1:filler,:)];
              for ji = 1:filler
              Routput{size(Routput,1)-ji+1,2} = '';
              Routput{size(Routput,1)-ji+1,3} = '';
              Routput{size(Routput,1)-ji+1,4} = '';
              end
            elseif size(Loutput,1) < size(Routput,1)
              filler = size(Routput,1) - size(Loutput,1);
              Loutput = [Loutput;Loutput(1:filler,:)];
              for ji = 1:filler
              Loutput{size(Loutput,1)-ji+1,2} = '';
              Loutput{size(Loutput,1)-ji+1,3} = '';
              Loutput{size(Loutput,1)-ji+1,4} = '';
              end
            end
            output = [output;Loutput,Routput];
        else
            
        % both ...
            index = ((IDX(:,1)==nex) + (IDX(:,2)==3)) == 2;
            output = data2string(output,measures(index,:),type);
            
        end
            
    end
    
    % write file   
        xlswrite([pwd,filesep,subject_name,'_ACLReport.xls'],output,'Sheet1','A1');
    
end

function output = data2string(output,data,type)
    
    if isempty(data)
        
        output{size(output,1)+1,2} = type;
        output{size(output,1)  ,3} = 'No files found';
        output{size(output,1)  ,4} = '';
        output{size(output,1)+1,2} = '';
        
    else

        output{size(output,1)+1,2} = type;
        p_trial = '';

        % overall summary
        for j = 1:size(data,1)

            if ~strcmp(p_trial,data.trial{j})
            output{size(output,1)+1,2} = data.trial{j};
                               p_trial = data.trial{j};
            else
            output{size(output,1)+1,2} = '';
            end

            output{size(output,1)  ,3} = data.measures{j};

                values = double(data(j,6:end));

            for n = 1:size(values,2)

                if isnan(values(n))
                output{size(output,1)  ,n+3} = '';    
                else
                output{size(output,1)  ,n+3} = num2str(values(n));
                end
                
                if j == 1
                    output{size(output,1)-1,n+3} = data.Properties.VarNames{n+5};
                end
            end
        end

        output{size(output,1)+1,2} = '';
        
        % mean summary
        if size(data,1) ~= max(data.VarIDX)
   
        output{size(output,1)+1,2} = 'mean';

        for j = 1:max(data.VarIDX)

            cdata = data(data.VarIDX == j,:);

            if j ~= 1
                output{size(output,1)+1,2} = '';
            end

            if size(cdata,1) == 1
                break
            else
                output{size(output,1),3} = cdata.measures{1};

                    values = double(cdata(:,6:end));

                for n = 1:size(values,2)
                    if sum(isnan(values(:,n))) == 0
                        output{size(output,1)  ,n+3} = num2str(mean(values(:,n)));   
                    else
                        val = values(~isnan(values(:,n)),n);
                        if isempty(val)
                            output{size(output,1)  ,n+3} = '';
                        else
                            output{size(output,1)  ,n+3} = num2str(mean(val));
                        end
                    end
                end            
            end
        end
        end
output{size(output,1)+1,2} = ''; %added by Kat 17.10.16 to add extra blank row after mean in spreadsheet after mean
    end
end

function sout = mergestruct(varargin)
%MERGESTRUCT Merge structures with unique fields.
%   Copyright 2009 The MathWorks, Inc.
% Start with collecting fieldnames, checking implicitly 
% that inputs are structures.
 fn = [];
    for k = 1:nargin
        try
            fn = [fn ; fieldnames(varargin{k})]; %#ok<AGROW>
        catch MEstruct
            throw(MEstruct)
        end
    end
    
    % Make sure the field names are unique.
    if length(fn) ~= length(unique(fn))
        error('mergestruct:FieldsNotUnique',...
            'Field names must be unique');
    end
    
    % Now concatenate the data from each struct.  Can't use 
    % structfun since input structs may not be scalar.
    c = [];
    for k = 1:nargin
        try
            c = [c ; struct2cell(varargin{k})]; %#ok<AGROW>
        catch MEdata
            throw(MEdata);
        end
    end
    
    % Construct the output.
    sout = cell2struct(c, fn, 1);

end

function BW = getBodyWeight(Modelfile)

   % open file
   fileID = fopen(Modelfile,'r');
   % read file line by line
   income = textscan(fileID,'%s','delimiter',',');
   income = income{1,1};
   % detect note section and the describtion
   fclose(fileID);
   
   for n = 1:size(income,1)
      if ~isempty(strfind(income{n},'Bodymass')) 
          BW = income{n};
          BW = str2double(BW(13:end));
          break
      end
   end

end