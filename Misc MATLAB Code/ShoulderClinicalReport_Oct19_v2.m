% Creates clinical report for shoulder lab testing
% Requires dataset to have ben reconstructed, labelled and modelled. 
% Looks for exercises called CMPU, BDL, PJ and JT that don't have 'exclude'
% in the Notes column

% Author: Kat Daniels
% Date: 07/10/19
clear
figs = 1;
 % Go to directory containing this function
%  origin = which('ShoulderClinicalReport_Oct19');
%  cd(origin(1:(max(strfind(origin,filesep)))-1))
 
 % Select Session folder
 origin = uigetdir('C:\', 'Select the Session folder');
 %origin = uigetdir('C:\Users\kdaniels\Desktop\ViconData\Performance', 'Select the Session folder');
 %origin = 'C:\Users\kdaniels\Desktop\ViconData\Performance\NOR 29041994\New Session';
 %origin = 'C:\Users\kdaniels\Desktop\ViconData\Performance\CT 24041988\New Session 1';
 %origin = 'E:\ProcessedStudyData\ShoulderNormative\LabData\AOC 22111987\New Session 1';
 %origin = 'C:\Users\kdaniels\Desktop\ViconData\Performance\LOC 204186 T1\New Session';
 %origin = 'C:\Users\kdaniels\Desktop\ViconData\Performance\MW 051119\New Session';
%origin = 'C:\Users\kdaniels\Desktop\ViconData\Performance\AD 151499 T1\New Session';

  
 % Make the user select the injured side. 1=left; 2= right
 [injside,~] = listdlg('ListString', {'Left', 'Right'},'SelectionMode','single', 'PromptString', 'Select the injured side',...
         'ListSize', [150,50], 'CancelString', 'Quit'); 
 
 % Get list of all c3d files in folder
 c3ds = dirPlus(origin, 'Depth', 0, 'FileFilter', '.c3d');
 
 % Get list of all enf files in folder
 enfs = dirPlus(origin, 'Depth', 0, 'FileFilter', '.enf');
 enftool = EnfToolbox;
 
 
 % ---------------------------------------------------------------------------------------------
 % ---------- Test if JPS present and only worry about Static Trial and modelling if so ---------
 % ---------------------------------------------------------------------------------------------
 
  % ---------------------------------------------------------------------  
 % ------ Get JPS results ---------------------------------------------- 
 % ---------------------------------------------------------------------
 JTc3ds = c3ds(contains(c3ds, {'JTLeft', 'JTRight'})); % Get JT (JPS open chain) trial c3ds
 
 if  ~isempty(JTc3ds) % If JT trials are present... 1==0;%

 % ---------------------------------------------------------------------  
 % ---------- Check Static Trial is modelled and get body mass ---------
 % ---------------------------------------------------------------------

 static = char(c3ds(contains(c3ds, 'Static')));
 if isempty(static)
     disp('STATIC TRIAL ISSUE: Cannot identify a static trial - is the trial name spelled correctly?')
     return
 end
 
 if size(static,1) ~=1 % if more than one static trial present, get user to select 
     staticcell = cellstr(static);
     [IDXstatic,tf] = listdlg('ListString', staticcell,'SelectionMode','single', 'PromptString', 'Multiple static trials identified - please select one',...
         'ListSize', [600,100], 'CancelString', 'Quit');
     if tf
         static = char(staticcell(IDXstatic));
     else
         return
     end
 end
 
 c3dfile = btkReadAcquisition(static); % read in c3d file
 try
    markers = btkGetMarkers(c3dfile);
    markers.RGH; % Test whether file is modelled by trying to get modelled joint centre data
    markers.LGH;
 catch
     disp('STATIC TRIAL ISSUE: Unable to get shoulder joint centre modelled markers - is the trial modelled?')
     btkCloseAcquisition(c3dfile)
     return
 end
 md = btkGetMetaData(c3dfile) ; % get meta data
 
 try
    bodymass = md.children.PROCESSING.children.Bodymass.info.values; % Get subject body mass (entered into model when setting up for data collection) 
 catch
    disp('STATIC TRIAL ISSUE: Unable to get body mass - has it been entered?')
    btkCloseAcquisition(c3dfile)
    return
 end
 btkCloseAcquisition(c3dfile)

 
 cleft=0; % initialise counter for JPS left figures
 cright=0; % initialise counter for JPS right figures
 figure(1); clf; % initialise figure to plot JPS left graphs
 figure(2); clf; % initialise figure to plot JPS right graphs
 subplotrows = 4  ; % set number of rows to plot JPS subplots in
 subplotcols = 3  ; % set number of columns to plot JPS subplots in
 
 matchError_left_high = []; % initialise variables
 matchError_left_low = [];
 matchError_right_high = [];
 matchError_right_low = [];
   
 matchError_left_high_names = []; % initialise variables
 matchError_left_low_names = [];
 matchError_right_high_names = [];
 matchError_right_low_names = [];
 
 
  for i = 1:length(JTc3ds) % for each JT trial
  trial = char(JTc3ds(i)); % get trial name
  trialminusfiletype = trial(1:end-4);
  filesepIDX = strfind(trialminusfiletype, filesep);
  trialname = trialminusfiletype(filesepIDX(end)+1:end) ; % get trial name without path - used for screen displays
  trialenf = char(enfs(contains(enfs,trialminusfiletype))); % get corresponding enf file  
  % get Note and Description from enf
        if isempty(trialenf)
            sprintf('JPS ISSUE: No enf file found for trial %s - trial will not be processed any further', trialname)
            continue
        else          
            [~,DESC] = enftool.GetNote(trialenf,'','DESCRIPTION=');
            DESC = DESC(~isspace(DESC)); % strip any whitespaces
            [~,NOTE] = enftool.GetNote(trialenf,'','NOTES=');
            NOTE = NOTE(~isspace(NOTE)); % strip any whitespaces
        end
      
  if strcmp(NOTE, 'exclude')
    sprintf('JPS trial %s has been excluded', trialname)  
    continue % skip if trial has been excluded
  end
 
  if ~strcmp(DESC, '=M') && ~strcmp(DESC, '=L') && ~strcmp(DESC, '=H')
    sprintf('JPS ISSUE: Nexus "Note" defining target height for %s is missing or contains a typo - trial will be skipped', trialname)
    continue
  end
 
   c3dfile = btkReadAcquisition(trial); % read in c3d file 
   markers = btkGetMarkers(c3dfile) ;% Get marker position data
   
    % Check which direction the patient is facing by relative y positions of
    % RASI and RPSI
    if nanmean(markers.STRN(:,2)) > nanmean(markers.T10(:,2))
        sprintf('JPS ISSUE: Patient is facing the wrong way in %s', trialname)
        continue
    end

    % Get side of interest (left/right). Identified from filename.
    if contains(trial, 'Right'); sidei= 'Right';
      elseif contains(trial, 'Left'); sidei = 'Left';
      else sprintf('JPS ISSUE: unable to identify whether %s is a Left or Right trial - check filename', trialname)
          continue
    end
     
%%%%%%%% Define markers and angles of relevance based on identified side %%
% 1) Get required model parameters and joint positions from Nexus:
try    
    if strcmp(sidei, 'Right') 
        cright = cright+1; % increment counter    
        if isfield(markers, 'RHUO') % If PiG model outputs
            HUO = markers.RHUO; % elbow joint coordinates
            CLO = markers.RCLO; % shoulder joint coordinates              
        elseif isfield(markers, 'REJC') ... else if new model outputs
            HUO = markers.REJC; % elbow joint coordinates
            CLO = markers.RGH; % shoulder joint coordinates 
        end
                
    elseif strcmp(sidei, 'Left') 
        cleft = cleft+1; % increment counter 
        if isfield(markers, 'LHUO') % If PiG model outputs
            HUO = markers.LHUO;
            CLO = markers.LCLO;  
        elseif isfield(markers, 'LEJC') ... else if new model outputs
            HUO = markers.LEJC; % elbow joint coordinates
            CLO = markers.LGH; % shoulder joint coordinates 
        end      
    end
   
catch
    sprintf('JPS ISSUE: No/incomplete modelled data found for trial %s', trialname)  
end
  
% 2) Calculate additional parameters:
% a) Humerus angle in global
HUMdelta = CLO - HUO; 
[sagtheta,rho] = cart2pol(HUMdelta(:,2), HUMdelta(:,3)) ;  
sagthetadeg = rad2deg(sagtheta); % convert to degrees
sagthetadeg2 = 90 - sagthetadeg; % sagittal plane angle (raising/extending arm = positive)
 
HUMglob_sag = sagthetadeg2; %  Angle of humerus in global coordinates in sagittal plane 

%%%%%%% Identify indices to position and reposition %%%
% 1) Identify sections to search for plateaus within:
Ydist = CLO(:,2) - HUO(:,2);
Zdist = CLO(:,3) - HUO(:,3);
[theta, rho] = cart2pol(Ydist, Zdist); % get elevation angle (theta) of arm from prox humerus prox to dist radius.

minang = min(theta); startang = max(theta); %find most-elevated and least-elevated arm angle for the trial 
threshProportion = 0.6 ; % set threshold to segment plateau region search area as a proportion of the way from start to max
getThresh = startang-(threshProportion*(startang - minang)) ; % get threshold value
sectionsIDX = find(theta<getThresh); % get indices to frames where theta is beyond the threshold
gapsize = diff(sectionsIDX);
split = find(gapsize>300); 
pos1_start = sectionsIDX(1);
pos1_end = sectionsIDX(split);
pos2_start = pos1_end + gapsize(split);
pos2_end = sectionsIDX(end);    
   
% 2) Search within sections for stable regions:
a = smooth(theta(pos1_start:pos1_end),20);
b = [nan; diff(a)];
c = find(b<0.0003 & b>-0.0003);

d = smooth(theta(pos2_start:pos2_end),20);
e = [nan; diff(d)];
f = find(e<0.0003 & e>-0.0003);

pos1_unedited = [c+pos1_start theta(c+pos1_start)]; %  idx, val
pos2_unedited = [f+pos2_start theta(f+pos2_start)]; %  idx, val
    
% 3) Remove any points outside xsd of mean for identified position and
% reposition:
x = 1.5; % Set number of SDs outside of which points are excluded
pos1mean = mean(pos1_unedited(:,2)); pos2mean = mean(pos2_unedited(:,2)); 
pos1sd = std(pos1_unedited(:,2)); pos2sd = std(pos2_unedited(:,2));

IDX_pos1withinrange = find(pos1_unedited(:,2)<(pos1mean+pos1sd*x) & pos1_unedited(:,2)>(pos1mean-pos1sd*x));
IDX_pos2withinrange = find(pos2_unedited(:,2)<(pos2mean+pos2sd*x) & pos2_unedited(:,2)>(pos2mean-pos2sd*x));

pos1_edited = pos1_unedited(IDX_pos1withinrange,:); pos2_edited = pos2_unedited(IDX_pos2withinrange,:); % idx, val (outliers removed)
     
HUMglob_sag_pos1mean =  mean(HUMglob_sag(pos1_edited(:,1))); % get arm angle in first position 
HUMglob_sag_pos2mean =  mean(HUMglob_sag(pos2_edited(:,1))); % get arm angle in reposition 
HUMglob_sag_ma = abs(HUMglob_sag_pos2mean - HUMglob_sag_pos1mean); % calculate absolute matching error as the absolute difference between them

if strcmp(sidei, 'Right') & strcmp(DESC, '=L')
    matchError_right_low = [matchError_right_low; HUMglob_sag_ma];
    matchError_right_low_names = [matchError_right_low_names; {trialname}];
elseif strcmp(sidei, 'Right') & strcmp(DESC, '=H')
    matchError_right_high = [matchError_right_high; HUMglob_sag_ma];
    matchError_right_high_names = [matchError_right_high_names; {trialname}];
elseif 'Left' & strcmp(DESC, '=L')
    matchError_left_low = [matchError_left_low; HUMglob_sag_ma];
    matchError_left_low_names = [matchError_left_low_names; {trialname}];
elseif 'Left' & strcmp(DESC, '=H')
    matchError_left_high = [matchError_left_high; HUMglob_sag_ma];
    matchError_left_high_names = [matchError_left_high_names; {trialname}];
end
 
if figs == 1
    if strcmp(sidei, 'Right') 
        figure(1);
        f = gcf;
        f.Units = 'centimeters';
        f.OuterPosition = [20 5 15 15];
        subplot(subplotrows,subplotcols,cright)
        hold on
        title(trialname)
        plot(HUMglob_sag) 
        scatter(sectionsIDX, HUMglob_sag(sectionsIDX))
        scatter([pos1_start pos1_end], [HUMglob_sag(pos1_start),HUMglob_sag(pos1_end)], 'g*')
        scatter([pos2_start pos2_end], [HUMglob_sag(pos2_start),HUMglob_sag(pos2_end)], 'c*')       
        scatter(pos1_unedited(:,1), HUMglob_sag(pos1_unedited(:,1)))
        scatter(pos2_unedited(:,1), HUMglob_sag(pos2_unedited(:,1)))    
        scatter(pos1_edited(:,1), HUMglob_sag(pos1_edited(:,1)), 'filled')
        scatter(pos2_edited(:,1), HUMglob_sag(pos2_edited(:,1)), 'filled') 
        hold off
    elseif strcmp(sidei, 'Left') 
        figure(2);
        f = gcf;
        f.Units = 'centimeters';
        f.OuterPosition = [3 5 15 15];
        subplot(subplotrows,subplotcols,cleft)
        hold on
        title(trialname)
        plot(HUMglob_sag) 
        scatter(sectionsIDX, HUMglob_sag(sectionsIDX))
        scatter([pos1_start pos1_end], [HUMglob_sag(pos1_start),HUMglob_sag(pos1_end)], 'g*')
        scatter([pos2_start pos2_end], [HUMglob_sag(pos2_start),HUMglob_sag(pos2_end)], 'c*')       
        scatter(pos1_unedited(:,1), HUMglob_sag(pos1_unedited(:,1)))
        scatter(pos2_unedited(:,1), HUMglob_sag(pos2_unedited(:,1)))    
        scatter(pos1_edited(:,1), HUMglob_sag(pos1_edited(:,1)), 'filled')
        scatter(pos2_edited(:,1), HUMglob_sag(pos2_edited(:,1)), 'filled') 
        hold off           
    end     
end 
btkCloseAcquisition(c3dfile) % close trial c3d file  
  end

  pause % wait for keypress
  % Make the user confirm the JPS plot results are correct before continuing.
 [allOK,~] = listdlg('ListString', {'Yes', 'No'},'SelectionMode','single', 'PromptString', 'All OK to continue?',...
  'ListSize', [150,50], 'CancelString', 'Quit');  

 close all % close the figures

 if allOK == 2 
     return % Quit function if not
 end

  % Calculate mean matching errors
  matchError_left_low_mean = nanmean(matchError_left_low);
  matchError_left_high_mean = nanmean(matchError_left_high);
  matchError_right_low_mean = nanmean(matchError_right_low);
  matchError_right_high_mean = nanmean(matchError_right_high);
  
else
  disp('JPS ISSUE: No JPS trials identified. If this is a surprise then go and check trial naming!')
  
  % Just get body mass from the first trial 
  firstTrial = char(c3ds(1));
  trialToUse = btkReadAcquisition(firstTrial); % read in c3d file
  md = btkGetMetaData(trialToUse) ; % get meta data  
  try
      bodymass = md.children.PROCESSING.children.Bodymass.info.values; % Get subject body mass (entered into model when setting up for data collection) 
      btkCloseAcquisition(trialToUse);
  catch
      disp('Unable to get body mass - has it been entered?')
  end

 end

 
 % ---------------------------------------------------------------------------------   
 % ------ Get CMPU results ---------------------------------------------- 
 % ---------------------------------------------------------------------------------  
  CMPUc3ds = c3ds(contains(c3ds, 'CMPU')); % Get CMPU trial c3ds
  if isempty(CMPUc3ds)
     disp('CMPU ISSUE: No CMPU trials identified. If this is a surprise then go and check trial naming!')
  else
 
CMPU_PEAK_TO_R = []; % initialise variables
CMPU_PEAK_TO_L = [];
CMPU_PEAK_LAN_R = [];
CMPU_PEAK_LAN_L = [];
CMPU_ECC2_IMP_R = [];
CMPU_ECC2_IMP_L = [];
CMPU_CONC_IMP_R = [];
CMPU_CONC_IMP_L = [];
CMPU_JH = [];
 
  for i = 1:length(CMPUc3ds) % for each CMPU trial
  trial = char(CMPUc3ds(i)); % get trial name
  trialminusfiletype = trial(1:end-4);
  filesepIDX = strfind(trialminusfiletype, filesep);
  trialname = trialminusfiletype(filesepIDX(end)+1:end) ; % get trial name without path - used for screen displays
  trialenf = char(enfs(contains(enfs,trialminusfiletype))); % get corresponding enf file  
  % get Note and Description from enf
        if isempty(trialenf)
            sprintf('CMPU ISSUE: No enf file found for trial %s - trial will not be processed any further', trialname)
            continue
        else          
            [~,DESC] = enftool.GetNote(trialenf,'','DESCRIPTION=');
            DESC = DESC(~isspace(DESC)); % strip any whitespaces
            [~,NOTE] = enftool.GetNote(trialenf,'','NOTES=');
            NOTE = NOTE(~isspace(NOTE)); % strip any whitespaces
        end      
  if strcmp(NOTE, 'exclude')
    sprintf('CMPU trial %s has been excluded', trialname)  
    continue % skip if trial has been excluded
  end
 
   c3dfile = btkReadAcquisition(trial); % read in c3d file 
   markers = btkGetMarkers(c3dfile) ;% Get marker position data
   analogs = btkGetAnalogs(c3dfile); % Get force data
   fpfreq = btkGetAnalogFrequency(c3dfile); % get forceplate frequency
   ratio = btkGetAnalogSampleNumberPerFrame(c3dfile); % ratio analog freq to modelled freq
   Fz1 = analogs.Force_Fz1*-1 ; 
   Fz2 = analogs.Force_Fz2*-1 ; 
   force = [Fz1 Fz2];   
   GRF = Fz1 + Fz2;
   GRF_Right = Fz1; 
   GRF_Left = Fz2;   

try    
   c7z = markers.C7(:,3); 
catch
    sprintf('CMPU ISSUE: Unable to get C7 marker for trial %s - check marker labelling', trialname)  
end
        takeoff   = find(GRF < 30,1,'first');
        [~,threshold] = max(GRF(takeoff:end));
        threshold = threshold + takeoff - 1;
        impact  = find(GRF(1:threshold) < 30,1,'last');
       
        %Get peak forces
        peak_to_R = max(GRF_Right(1:takeoff));
        peak_to_L = max(GRF_Left(1:takeoff));        
        peak_lan_R = max(GRF_Right(impact:end));
        peak_lan_L = max(GRF_Left(impact:end));        
        
        % Get jump height from C7 marker displacement
        c7peakflight = max(c7z(floor(takeoff/ratio):floor(impact/ratio)));
        c7to = c7z(floor(takeoff/ratio));
        jumph = (c7peakflight - c7to) ./10; % in cm
        
        try % put all the impulse calculations in a try-catch conditional as not needed for the normal clinical report
        % Get impulses
            % Get end of take-off eccentric phase based on C7 marker position: 
            c7z_to = c7z(1:floor(takeoff/ratio)); %... c7 z position from start to take-off  
            ecc_end = [nan;diff(c7z_to)]; ecc_end = ecc_end > 0; ecc_end = [nan;diff(ecc_end)]; % get - to + turning points
            IDXecc_end = find(ecc_end==1, 1, 'last').*ratio ; % analogue frame index to end of eccentric phase in take-off
        
           % Get start of eccentric phase 2 (braking phase) based on velocity of neck marker:
            c7z_to_vv = [nan; diff(c7z_to)]; % diff to get vertical velocity  
            ecc2_start = [nan;diff(c7z_to_vv)]; ecc2_start = ecc2_start > 0; ecc2_start = [nan;diff(ecc2_start)]; % get - to + turning points
            ecc2_start(c7z_to_vv>0) = 0; % Exclude any turning points where velocity is not negative (downwards)
            IDXecc2_start = find(ecc2_start==1, 1, 'last').*ratio; % analogue frame index to start of eccentric phase 2 in take-off 
         
            ecc2_force_FP1 = Fz1(IDXecc2_start:IDXecc_end) ; 
            ecc2_impF1 = 1./fpfreq * (trapz(ecc2_force_FP1)); % get eccentric phase 2 impulse FP1 
            ecc2_force_FP2 = Fz2(IDXecc2_start:IDXecc_end) ; 
            ecc2_impF2 = 1./fpfreq * (trapz(ecc2_force_FP2)); % get eccentric phase 2 impulse FP2    
            
          % Get take-off concentric phase impulses    
            conc_force_FP1 = Fz1(IDXecc_end:takeoff) ; 
            conc_impF1 = 1./fpfreq * (trapz(conc_force_FP1)); % get take-off concentric phase impulse FP1 
            conc_force_FP2 = Fz2(IDXecc_end:takeoff) ; 
            conc_impF2 = 1./fpfreq * (trapz(conc_force_FP2)); % get take-off concentric phase impulse FP2    
         
           % Express impulses relative to body mass  
            ecc2_imp_R = ecc2_impF1/bodymass;
            ecc2_imp_L = ecc2_impF2/bodymass;
            conc_imp_R = conc_impF1/bodymass;
            conc_imp_L = conc_impF2/bodymass;            
            
        catch
            sprintf('CMPU ISSUE: Issue with impulse calculations for trial %s - check if needed. All returned as NaN.', trialname)  
            ecc2_imp_R = nan;
            ecc2_imp_L = nan;
            conc_imp_R = nan;
            conc_imp_L = nan;
        end
        
         % Express peak forces relative to body mass  
         peak_to_R = peak_to_R/bodymass;
         peak_to_L = peak_to_L/bodymass;
         peak_lan_R = peak_lan_R/bodymass;
         peak_lan_L = peak_lan_L/bodymass;
         
    CMPU_PEAK_TO_R = [CMPU_PEAK_TO_R; peak_to_R]; 
    CMPU_PEAK_TO_L = [CMPU_PEAK_TO_L; peak_to_L];
    CMPU_PEAK_LAN_R = [CMPU_PEAK_LAN_R; peak_lan_R];
    CMPU_PEAK_LAN_L = [CMPU_PEAK_LAN_L; peak_lan_L];
    CMPU_ECC2_IMP_R = [CMPU_ECC2_IMP_R; ecc2_imp_R];
    CMPU_ECC2_IMP_L = [CMPU_ECC2_IMP_L; ecc2_imp_L];
    CMPU_CONC_IMP_R = [CMPU_CONC_IMP_R; conc_imp_R];
    CMPU_CONC_IMP_L = [CMPU_CONC_IMP_L; conc_imp_L];    
    CMPU_JH = [CMPU_JH; jumph];
      
  end
  % Calculate mean values
  CMPU_PEAK_TO_R_mean = nanmean(CMPU_PEAK_TO_R);
  CMPU_PEAK_TO_L_mean = nanmean(CMPU_PEAK_TO_L);
  CMPU_PEAK_LAN_R_mean = nanmean(CMPU_PEAK_LAN_R);
  CMPU_PEAK_LAN_L_mean = nanmean(CMPU_PEAK_LAN_L);
  CMPU_ECC2_IMP_R_mean = nanmean(CMPU_ECC2_IMP_R);
  CMPU_ECC2_IMP_L_mean = nanmean(CMPU_ECC2_IMP_L);
  CMPU_CONC_IMP_R_mean = nanmean(CMPU_CONC_IMP_R);
  CMPU_CONC_IMP_L_mean = nanmean(CMPU_CONC_IMP_L);
  CMPU_JH_mean = nanmean(CMPU_JH);
  
  % Identify best trial
  [~,IDX_maxCMPU] = max(CMPU_JH);
  CMPUbesttrial = CMPUc3ds(IDX_maxCMPU);
  
  CMPUbesttrial = char(CMPUbesttrial); % get trial name
  CMPUbesttrial = CMPUbesttrial(1:end-4);
  filesepIDX1 = strfind(CMPUbesttrial, filesep);
  CMPUbesttrial = CMPUbesttrial(filesepIDX1(end)+1:end)  ;
  
  CMPU_PEAK_TO_R_best = CMPU_PEAK_TO_R(IDX_maxCMPU);
  CMPU_PEAK_TO_L_best = CMPU_PEAK_TO_L(IDX_maxCMPU);
  CMPU_PEAK_LAN_R_best = CMPU_PEAK_LAN_R(IDX_maxCMPU);
  CMPU_PEAK_LAN_L_best = CMPU_PEAK_LAN_L(IDX_maxCMPU);
  CMPU_ECC2_IMP_R_best = CMPU_ECC2_IMP_R(IDX_maxCMPU);
  CMPU_ECC2_IMP_L_best = CMPU_ECC2_IMP_L(IDX_maxCMPU);
  CMPU_CONC_IMP_R_best = CMPU_CONC_IMP_R(IDX_maxCMPU);
  CMPU_CONC_IMP_L_best = CMPU_CONC_IMP_L(IDX_maxCMPU);
  CMPU_JH_best = CMPU_JH(IDX_maxCMPU);
    
  end

% ---------------------------------------------------------------------------------   
 % ------ Get BDL results ---------------------------------------------- 
 % ---------------------------------------------------------------------------------  
  BDLc3ds = c3ds(contains(c3ds, 'BDL')); % Get CMPU trial c3ds
  if isempty(BDLc3ds)
     disp('BDL ISSUE: No BDL trials identified. If this is a surprise then go and check trial naming!')
  else
  

BDL_PEAK_LAN_R = []; % initialise variables
BDL_PEAK_LAN_L = [];
 
  for i = 1:length(BDLc3ds) % for each CMPU trial
  trial = char(BDLc3ds(i)); % get trial name
  trialminusfiletype = trial(1:end-4);
  filesepIDX = strfind(trialminusfiletype, filesep);
  trialname = trialminusfiletype(filesepIDX(end)+1:end) ; % get trial name without path - used for screen displays
  trialenf = char(enfs(contains(enfs,trialminusfiletype))); % get corresponding enf file  
  % get Note and Description from enf
        if isempty(trialenf)
            sprintf('BDL ISSUE: No enf file found for trial %s - trial will not be processed any further', trialname)
            continue
        else          
            [~,DESC] = enftool.GetNote(trialenf,'','DESCRIPTION=');
            DESC = DESC(~isspace(DESC)); % strip any whitespaces
            [~,NOTE] = enftool.GetNote(trialenf,'','NOTES=');
            NOTE = NOTE(~isspace(NOTE)); % strip any whitespaces
        end      
  if strcmp(NOTE, 'exclude')
    sprintf('BDL trial %s has been excluded', trialname)  
    continue % skip if trial has been excluded
  end
 
   c3dfile = btkReadAcquisition(trial); % read in c3d file 
   markers = btkGetMarkers(c3dfile) ;% Get marker position data
   analogs = btkGetAnalogs(c3dfile); % Get force data
   ratio = btkGetAnalogSampleNumberPerFrame(c3dfile); % ratio analog freq to modelled freq
   Fz1 = analogs.Force_Fz1*-1 ; 
   Fz2 = analogs.Force_Fz2*-1 ; 
   force = [Fz1 Fz2];   
   GRF = Fz1 + Fz2;
   GRF_Right = Fz1; 
   GRF_Left = Fz2;   

       impact  = find(GRF < 30,1,'first');
       
        %Get peak forces      
        peak_lan_R = max(GRF_Right(impact:end));
        peak_lan_L = max(GRF_Left(impact:end));        
              
        % Express peak forces relative to body mass  
        peak_lan_R = peak_lan_R/bodymass;
        peak_lan_L = peak_lan_L/bodymass;
         
    BDL_PEAK_LAN_R = [BDL_PEAK_LAN_R; peak_lan_R];
    BDL_PEAK_LAN_L = [BDL_PEAK_LAN_L; peak_lan_L];

  end
  % Calculate mean values
  BDL_PEAK_LAN_R_mean = nanmean(BDL_PEAK_LAN_R);
  BDL_PEAK_LAN_L_mean = nanmean(BDL_PEAK_LAN_L);

  
   % Identify best trial based on max sum of peak forces (L+R)
  [~,IDX_maxBDL] = max(sum([BDL_PEAK_LAN_R BDL_PEAK_LAN_L], 2));
  BDLbesttrial = BDLc3ds(IDX_maxBDL);
  
  BDLbesttrial = char(BDLbesttrial); % get trial name
  BDLbesttrial = BDLbesttrial(1:end-4);
  filesepIDX1 = strfind(BDLbesttrial, filesep);
  BDLbesttrial = BDLbesttrial(filesepIDX1(end)+1:end)  ;
   
  BDL_PEAK_LAN_R_best = BDL_PEAK_LAN_R(IDX_maxBDL);
  BDL_PEAK_LAN_L_best = BDL_PEAK_LAN_L(IDX_maxBDL);
   
  end
  
  
  
  % ---------------------------------------------------------------------------------   
 % ------ Get PressJump results ---------------------------------------------- 
 % ---------------------------------------------------------------------------------  
  PJc3ds = c3ds(contains(c3ds, 'PJ')); % Get CMPU trial c3ds
  if isempty(PJc3ds)
     disp('PJ ISSUE: No PJ trials identified. If this is a surprise then go and check trial naming!')
  else
 
PJ_PEAK_TO_R = []; % initialise variables
PJ_PEAK_TO_L = [];
PJ_CONC_IMP_R = [];
PJ_CONC_IMP_L = [];
PJ_JH = [];
 
  for i = 1:length(PJc3ds) % for each PJ trial
  trial = char(PJc3ds(i)); % get trial name
  trialminusfiletype = trial(1:end-4);
  filesepIDX = strfind(trialminusfiletype, filesep);
  trialname = trialminusfiletype(filesepIDX(end)+1:end) ; % get trial name without path - used for screen displays
  trialenf = char(enfs(contains(enfs,trialminusfiletype))); % get corresponding enf file  
  % get Note and Description from enf
        if isempty(trialenf)
            sprintf('PJ ISSUE: No enf file found for trial %s - trial will not be processed any further', trialname)
            continue
        else          
            [~,DESC] = enftool.GetNote(trialenf,'','DESCRIPTION=');
            DESC = DESC(~isspace(DESC)); % strip any whitespaces
            [~,NOTE] = enftool.GetNote(trialenf,'','NOTES=');
            NOTE = NOTE(~isspace(NOTE)); % strip any whitespaces
        end      
  if strcmp(NOTE, 'exclude')
    sprintf('PJ trial %s has been excluded', trialname)  
    continue % skip if trial has been excluded
  end
 
   c3dfile = btkReadAcquisition(trial); % read in c3d file 
   markers = btkGetMarkers(c3dfile) ;% Get marker position data
   analogs = btkGetAnalogs(c3dfile); % Get force data
   ratio = btkGetAnalogSampleNumberPerFrame(c3dfile); % ratio analog freq to modelled freq
   Fz1 = analogs.Force_Fz1*-1 ; 
   Fz2 = analogs.Force_Fz2*-1 ; 
   force = [Fz1 Fz2];   
   GRF = Fz1 + Fz2;
   GRF_Right = Fz1; 
   GRF_Left = Fz2;   

try    
   c7z = markers.C7(:,3); 
catch
    sprintf('PJ ISSUE: Unable to get C7 marker for trial %s - check marker labelling', trialname)  
end
        takeoff   = find(GRF < 30,1,'first');
        [~,threshold] = max(GRF(takeoff:end));
        threshold = threshold + takeoff - 1;
        impact  = find(GRF(1:threshold) < 30,1,'last');
       
        %Get peak forces
        peak_to_R = max(GRF_Right(1:takeoff));
        peak_to_L = max(GRF_Left(1:takeoff));        
        peak_lan_R = max(GRF_Right(impact:end));
        peak_lan_L = max(GRF_Left(impact:end));        
        
        % Get jump height from C7 marker displacement
        c7peakflight = max(c7z(floor(takeoff/ratio):floor(impact/ratio)));
        c7to = c7z(floor(takeoff/ratio));
        jumph = (c7peakflight - c7to) ./10; % in cm
        
        try % put all the impulse calculations in a try-catch conditional as not needed for the normal clinical report
        % Get impulses
            % Get end of take-off eccentric phase based on C7 marker position: 
            c7z_to = c7z(1:floor(takeoff/ratio)); %... c7 z position from start to take-off  
            ecc_end = [nan;diff(c7z_to)]; ecc_end = ecc_end > 0; ecc_end = [nan;diff(ecc_end)]; % get - to + turning points
            IDXecc_end = find(ecc_end==1, 1, 'last').*ratio ; % analogue frame index to end of eccentric phase in take-off
                   
          % Get take-off concentric phase impulses    
            conc_force_FP1 = Fz1(IDXecc_end:takeoff) ; 
            conc_impF1 = 1./fpfreq * (trapz(conc_force_FP1)); % get take-off concentric phase impulse FP1 
            conc_force_FP2 = Fz2(IDXecc_end:takeoff) ; 
            conc_impF2 = 1./fpfreq * (trapz(conc_force_FP2)); % get take-off concentric phase impulse FP2    
         
           % Express impulses relative to body mass  
            conc_imp_R = conc_impF1/bodymass;
            conc_imp_L = conc_impF2/bodymass;            
            
        catch
            sprintf('PJ ISSUE: Issue with impulse calculations for trial %s - check if needed. All returned as NaN.', trialname)  
            conc_imp_R = nan;
            conc_imp_L = nan;
        end
        
         % Express peak forces relative to body mass  
         peak_to_R = peak_to_R/bodymass;
         peak_to_L = peak_to_L/bodymass;
         
    PJ_PEAK_TO_R = [PJ_PEAK_TO_R; peak_to_R]; 
    PJ_PEAK_TO_L = [PJ_PEAK_TO_L; peak_to_L];
    PJ_CONC_IMP_R = [PJ_CONC_IMP_R; conc_imp_R];
    PJ_CONC_IMP_L = [PJ_CONC_IMP_L; conc_imp_L];       
    PJ_JH = [PJ_JH; jumph]; 
  end
  % Calculate mean values
  PJ_PEAK_TO_R_mean = nanmean(CMPU_PEAK_TO_R);
  PJ_PEAK_TO_L_mean = nanmean(CMPU_PEAK_TO_L);
  PJ_CONC_IMP_R_mean = nanmean(CMPU_CONC_IMP_R);
  PJ_CONC_IMP_L_mean = nanmean(CMPU_CONC_IMP_L);
  CMPU_JH_mean = nanmean(CMPU_JH); 
  
  % Identify best trial
  [~,IDX_maxPJ] = max(PJ_JH);
  PJbesttrial = PJc3ds(IDX_maxPJ);
  
  PJbesttrial = char(PJbesttrial); % get trial name
  PJbesttrial = PJbesttrial(1:end-4);
  filesepIDX1 = strfind(PJbesttrial, filesep);
  PJbesttrial = PJbesttrial(filesepIDX1(end)+1:end)  ;
  
  PJ_PEAK_TO_R_best = PJ_PEAK_TO_R(IDX_maxPJ);
  PJ_PEAK_TO_L_best = PJ_PEAK_TO_L(IDX_maxPJ);
  PJ_CONC_IMP_R_best = PJ_CONC_IMP_R(IDX_maxPJ);
  PJ_CONC_IMP_L_best = PJ_CONC_IMP_L(IDX_maxPJ);
  PJ_JH_best = PJ_JH(IDX_maxPJ);
  end
  
  
   
 % ---------------------------------------------------------------------------------   
 % ------ Produce clinical report --------------------------------------------------
 % ---------------------------------------------------------------------------------  
  filler = {'','','';'','',''};
   %------------------------------------------------------------------------
 % Sheet 1: Standard clinical report 
  %-------------------------------------------------------------------------
 try
 CMPUlabels_Row1 = {'CMPU'; ''; 'Peak take-off force (N/kg)'; ''; ''; ''; 'Peak landing force (N/kg)'; ''; ''; ''; 'Jump height'};
 CMPUlabels_Row2 = {CMPUbesttrial; ''; 'Left'; 'Right'; '% Diff'; ''; 'Left'; 'Right'; '% Diff'; ''; 'cm' }; 
 if injside ==1 % if left side inj
     diff_TO = (CMPU_PEAK_TO_R_best - CMPU_PEAK_TO_L_best) ./CMPU_PEAK_TO_R_best .*100;
     diff_LAN = (CMPU_PEAK_LAN_R_best - CMPU_PEAK_LAN_L_best) ./CMPU_PEAK_LAN_R_best .*100;
 else % if right side inj
     diff_TO = (CMPU_PEAK_TO_L_best - CMPU_PEAK_TO_R_best) ./CMPU_PEAK_TO_L_best .*100;
     diff_LAN = (CMPU_PEAK_LAN_L_best - CMPU_PEAK_LAN_R_best) ./CMPU_PEAK_LAN_L_best .*100;
 end 
 CMPUdata = {''; ''; round(CMPU_PEAK_TO_L_best,1); round(CMPU_PEAK_TO_R_best,1); round(diff_TO,1); ''; round(CMPU_PEAK_LAN_L_best,1); round(CMPU_PEAK_LAN_R_best,1); round(diff_LAN,1); ''; round(CMPU_JH_best,1)};  
 CMPU_summary_outputs = [CMPUlabels_Row1 CMPUlabels_Row2 CMPUdata]; % CMPU results from best trial
 catch
 CMPU_summary_outputs ={'', '', ''};
 end
 
 try
 BDLlabels_Row1 = {'BDL'; ''; 'Peak landing force (N/kg)'; ''; ''; };
 BDLlabels_Row2 = {BDLbesttrial; ''; 'Left'; 'Right'; '% Diff';}; 
 if injside ==1 % if left side inj
     diff_LAN = (BDL_PEAK_LAN_R_best - BDL_PEAK_LAN_L_best) ./BDL_PEAK_LAN_R_best .*100;
 else % if right side inj
     diff_LAN = (BDL_PEAK_LAN_L_best - BDL_PEAK_LAN_R_best) ./BDL_PEAK_LAN_L_best .*100;
 end 
 BDLdata = {''; ''; round(BDL_PEAK_LAN_L_best,1); round(BDL_PEAK_LAN_R_best,1); round(diff_LAN,1);};  
 BDL_summary_outputs = [BDLlabels_Row1 BDLlabels_Row2 BDLdata]; % BDL results from best trial
 catch
 BDL_summary_outputs ={'', '', ''};
 end
 
 try
 PJlabels_Row1 = {'PJ'; ''; 'Peak take-off force (N/kg)'; ''; ''; ''; 'Jump height'};
 PJlabels_Row2 = {PJbesttrial; ''; 'Left'; 'Right'; '% Diff'; ''; 'cm' }; 
 if injside ==1 % if left side inj
     diff_TO = (PJ_PEAK_TO_R_best - PJ_PEAK_TO_L_best) ./PJ_PEAK_TO_R_best .*100;
 else % if right side inj
     diff_TO = (PJ_PEAK_TO_L_best - PJ_PEAK_TO_R_best) ./PJ_PEAK_TO_L_best .*100;
 end 
 PJdata = {''; ''; round(PJ_PEAK_TO_L_best,1); round(PJ_PEAK_TO_R_best,1); round(diff_TO,1); ''; round(PJ_JH_best,1)};  
 PJ_summary_outputs = [PJlabels_Row1 PJlabels_Row2 PJdata]; % CMPU results from best trial
 catch
 PJ_summary_outputs ={'', '', ''};
 end
 
 try
 JPSlabels_Row1 = {'JPS Matching Error (degrees)'; ''; 'Low'; ''; ''; ''; 'High'; ''; ''};
 JPSlabels_Row2 = {''; ''; 'Left'; 'Right'; 'Diff'; ''; 'Left'; 'Right'; 'Diff'}; 
 if injside ==1 % if left side inj
     matchError_low_diff = matchError_left_low_mean - matchError_right_low_mean;
     matchError_high_diff = matchError_left_high_mean - matchError_right_high_mean;
 else % if right side inj
     matchError_low_diff = matchError_right_low_mean - matchError_left_low_mean;
     matchError_high_diff = matchError_right_high_mean - matchError_left_high_mean;
 end  
 JPSdata = {''; ''; round(matchError_left_low_mean,1); round(matchError_right_low_mean,1); round(matchError_low_diff,1); ''; round(matchError_left_high_mean,1); round(matchError_right_high_mean,1); round(matchError_high_diff,1)};  
 JPS_summary_outputs = [JPSlabels_Row1 JPSlabels_Row2 JPSdata]; % JPS results from best trial
 catch
 JPS_summary_outputs = {'', '', ''}; 
 end

 SUMMARY_OUTPUTS = [CMPU_summary_outputs; filler; BDL_summary_outputs; filler; PJ_summary_outputs; filler; JPS_summary_outputs]; % Concatenate summary outputs for all exercises present
     
 %-------------------------------------------------------------------------
 % Sheet 2: Extended clinical report (+ impulses)
 %-------------------------------------------------------------------------
 try
 CMPUlabels_Row1 = {'CMPU'; ''; 'Peak take-off force (N/kg)'; ''; ''; ''; 'Eccentric decel impulse (N.s/kg)'; ''; ''; ''; 'Concentric impulse (N.s/kg)'; ''; ''; ''; 'Peak landing force (N/kg'; ''; ''; ''; 'Jump height'};
 CMPUlabels_Row2 = {CMPUbesttrial; ''; 'Left'; 'Right'; '% Diff'; ''; 'Left'; 'Right'; '% Diff'; ''; 'Left'; 'Right'; '% Diff'; '';'Left'; 'Right'; '% Diff'; ''; 'cm' }; 
 if injside ==1 % if left side inj
     diff_TOpf = (CMPU_PEAK_TO_R_best - CMPU_PEAK_TO_L_best) ./CMPU_PEAK_TO_R_best .*100;
     diff_TOeccimp = (CMPU_ECC2_IMP_R_best - CMPU_ECC2_IMP_L_best) ./CMPU_ECC2_IMP_R_best .*100;
     diff_TOconcimp = (CMPU_CONC_IMP_R_best - CMPU_CONC_IMP_L_best) ./CMPU_CONC_IMP_R_best .*100;
     diff_LANpf = (CMPU_PEAK_LAN_R_best - CMPU_PEAK_LAN_L_best) ./CMPU_PEAK_LAN_R_best .*100;
 else % if right side inj
     diff_TOpf = (CMPU_PEAK_TO_L_best - CMPU_PEAK_TO_R_best) ./CMPU_PEAK_TO_L_best .*100;
     diff_TOeccimp = (CMPU_ECC2_IMP_L_best - CMPU_ECC2_IMP_R_best) ./CMPU_ECC2_IMP_L_best .*100;
     diff_TOconcimp = (CMPU_CONC_IMP_L_best - CMPU_CONC_IMP_R_best) ./CMPU_CONC_IMP_L_best .*100;
     diff_LANpf = (CMPU_PEAK_LAN_L_best - CMPU_PEAK_LAN_R_best) ./CMPU_PEAK_LAN_L_best .*100;
 end 
 CMPUdata = {''; ''; round(CMPU_PEAK_TO_L_best,1); round(CMPU_PEAK_TO_R_best,1); round(diff_TOpf,1); ''; ...
     round(CMPU_ECC2_IMP_L_best,1); round(CMPU_ECC2_IMP_R_best,1); round(diff_TOeccimp,1); ''; ...
     round(CMPU_CONC_IMP_L_best,1); round(CMPU_CONC_IMP_R_best,1); round(diff_TOconcimp,1); ''; ...
     round(CMPU_PEAK_LAN_L_best,1); round(CMPU_PEAK_LAN_R_best,1); round(diff_LANpf,1); ''; round(CMPU_JH_best,1)};  
 CMPU_extended_outputs = [CMPUlabels_Row1 CMPUlabels_Row2 CMPUdata]; % CMPU results from best trial
 catch
 CMPU_extended_outputs ={'', '', ''};
 end
  
 
 try
 PJlabels_Row1 = {'PJ'; ''; 'Peak take-off force (N/kg)'; ''; ''; ''; 'Concentric impulse (N.s/kg)'; ''; ''; ''; 'Jump height'};
 PJlabels_Row2 = {PJbesttrial; ''; 'Left'; 'Right'; '% Diff'; '';'Left'; 'Right'; '% Diff'; ''; 'cm' }; 
 if injside ==1 % if left side inj
     diff_TOpf = (PJ_PEAK_TO_R_best - PJ_PEAK_TO_L_best) ./PJ_PEAK_TO_R_best .*100;
     diff_TOconcimp = (PJ_CONC_IMP_R_best - PJ_CONC_IMP_L_best) ./PJ_CONC_IMP_R_best .*100;
 else % if right side inj
     diff_TOpf = (PJ_PEAK_TO_L_best - PJ_PEAK_TO_R_best) ./PJ_PEAK_TO_L_best .*100;
     diff_TOconcimp = (PJ_CONC_IMP_L_best - PJ_CONC_IMP_R_best) ./PJ_CONC_IMP_L_best .*100;
 end 
 PJdata = {''; ''; round(PJ_PEAK_TO_L_best,1); round(PJ_PEAK_TO_R_best,1); round(diff_TOpf,1); '';...
     round(PJ_CONC_IMP_L_best,1); round(PJ_CONC_IMP_R_best,1); round(diff_TOconcimp,1); ''; round(PJ_JH_best,1)};  
 PJ_extended_outputs = [PJlabels_Row1 PJlabels_Row2 PJdata]; % CMPU results from best trial
 catch
 PJ_extended_outputs ={'', '', ''};
 end
 
 EXTENDED_OUTPUTS = [CMPU_extended_outputs; filler; BDL_summary_outputs; filler; PJ_extended_outputs; filler; JPS_summary_outputs]; % Concatenate summary outputs for all exercises present 
 
  %-------------------------------------------------------------------------   
 % Sheet 3: Full clinical report (all variables for all trials)
  %-------------------------------------------------------------------------

% JPS 
try
JPScellarray = {};
JPSlefthigh = [matchError_left_high_names num2cell(matchError_left_high)] ;
JPSrighthigh = [matchError_right_high_names num2cell(matchError_right_high)] ;
JPSleftlow = [matchError_left_low_names num2cell(matchError_left_low)] ;
JPSrightlow = [matchError_right_low_names num2cell(matchError_right_low)] ;

JPScellarray(1:size(JPSlefthigh,1),1:2) = JPSlefthigh;
JPScellarray(1:size(JPSrighthigh,1),4:5) = JPSrighthigh;
JPScellarray(1:size(JPSleftlow,1),7:8) = JPSleftlow;
JPScellarray(1:size(JPSrightlow,1),10:11) = JPSrightlow;

toprow = {'JPS Matching Error (degrees)', '','','','','','','','','',''};
secondrow = {'High Left','','','High Right', '','','Low Left', '','','Low Right','',};

JPSallResults = [toprow; secondrow; JPScellarray] ; % All JPS results in cell array ready to go
catch
   JPSallResults = {'','','','','','','','','','','',''}; 
end

% CMPU
try
CMPU_allresults = [CMPU_JH CMPU_PEAK_TO_L CMPU_ECC2_IMP_L CMPU_CONC_IMP_L CMPU_PEAK_LAN_L CMPU_PEAK_TO_R CMPU_ECC2_IMP_R CMPU_CONC_IMP_R CMPU_PEAK_LAN_R] ;
CMPU_allresults = num2cell(CMPU_allresults);
toprow = {'CMPU', '','','','','','','',''};
secondrow = {'JumpHeight (cm)' 'TakeOffPeakForce_Left (N/kg)' 'TakeOffEccDecImpulse_Left (N.s/kg)' 'TakeOffConcImpulse_Left (N.s/kg)' 'LandingPeakForce_Left (N/kg)'...
    'TakeOffPeakForce_Right (N/kg)' 'TakeOffEccDecImpulse_Right (N.s/kg)' 'TakeOffConcImpulse_Right (N.s/kg)' 'LandingPeakForce_Right (N/kg)'};
CMPUallResults = [toprow; secondrow; CMPU_allresults] ;% All CMPU results in cell array ready to go
catch
CMPUallResults = {'','','','','','','','',''};       
end

% BDL
try
BDL_allresults = [BDL_PEAK_LAN_L BDL_PEAK_LAN_R] ;
BDL_allresults = num2cell(BDL_allresults);
toprow = {'BDL', ''};
secondrow = {'LandingPeakForce_Left (N/kg)' 'LandingPeakForce_Right (N/kg)'};
BDLallResults = [toprow; secondrow; BDL_allresults] ;% All BDL results in cell array ready to go
catch
BDLallResults = {'',''};    
end

% PressJump
try
PJ_allresults = [PJ_JH PJ_PEAK_TO_L PJ_CONC_IMP_L PJ_PEAK_TO_R PJ_CONC_IMP_R] ;
PJ_allresults = num2cell(PJ_allresults);
toprow = {'PJ', '','','',''};
secondrow = {'JumpHeight (cm)' 'TakeOffPeakForce_Left (N/kg)' 'TakeOffConcImpulse_Left (N.s/kg)'...
    'TakeOffPeakForce_Right (N/kg)' 'TakeOffConcImpulse_Right (N.s/kg)' };
PJallResults = [toprow; secondrow; PJ_allresults] ;% All CMPU results in cell array ready to go
catch
PJallResults = {'','','','',''};    
end

% --> Combine
AllResultsCellArray = {};
AllResultsCellArray(1:size(CMPUallResults,1),1:size(CMPUallResults,2)) = CMPUallResults;
    currentnumrows = size(AllResultsCellArray,1)+2; 
AllResultsCellArray(currentnumrows:currentnumrows+size(BDLallResults,1)-1, 1:size(BDLallResults,2)) = BDLallResults;
    currentnumrows = size(AllResultsCellArray,1)+2; 
AllResultsCellArray(currentnumrows:currentnumrows+size(PJallResults,1)-1, 1:size(PJallResults,2)) = PJallResults;    
    currentnumrows = size(AllResultsCellArray,1)+2; 
AllResultsCellArray(currentnumrows:currentnumrows+size(JPSallResults,1)-1, 1:size(JPSallResults,2)) = JPSallResults;      
    
FULL_OUTPUTS =  AllResultsCellArray;

 %-------------------------------------------------------------------------   
 % Write to Excel
 %-------------------------------------------------------------------------
 filesepIDX = strfind(origin, filesep);
 patientname = origin(filesepIDX(end-1)+1:filesepIDX(end)-1) ; % get trial name without path - used for screen displays

 reportfilename = ['ShoulderLabReport_' patientname '_' date '.xlsx'];
 xlswrite(reportfilename,FULL_OUTPUTS,'FullOutputs')
 xlswrite(reportfilename,EXTENDED_OUTPUTS,'ExtendedOutputs')
 xlswrite(reportfilename,SUMMARY_OUTPUTS,'SummaryOutputs')

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
  