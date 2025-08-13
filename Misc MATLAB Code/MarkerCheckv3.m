% MarkerCheck2 written by K Daniels 01/02/2017 (beta)

% Script designed to be run as a Vicon Nexus 2 pipeline (tested in version
% 2.3) to write peak knee rotation and varus/valgus values to the trial
% history file for marker placement screening before commencing data capture 
% of test exercises.
% Requires marker T10 to be labelled and file to be modelled
% (PiG).

% v3 08/08/17


clear
disp('Starting checks...')
vicon = ViconNexus() ;

Subjects = vicon.GetSubjectNames() ; % Returns name of subject
Subjects = char(Subjects);

%--------------------
if 1==0 % section of code inserted in case needed for step detection in future
% get Toe marker trajectories
[LTOEx, LTOEy, LTOEz, LTOEe] = GetTrajectory( vicon, Subjects, 'LTOE' ) ;  % GetTrajectory get all frames of data for the trial for the specified marker. 
[RTOEx, RTOEy, RTOEz, RTOEe] = GetTrajectory( vicon, Subjects, 'RTOE' ) ;  % GetTrajectory get all frames of data for the trial for the specified marker. 
            % Returns
            %    x        = numerical(double) array, x-coordinates of the trajectory
            %    y        = numerical(double) array, y-coordinates of the trajectory
            %    z        = numerical(double) array, z-coordinates of the trajectory
            %    e        = logical array, T/F indication as to whether the data exists for each frame
end
%--------------------

% Identify index to frames where patient is inside the selected test area (marker
% T10)
try
[T10x, T10y, T10z, T10e] = GetTrajectory( vicon, Subjects, 'T10' ) ;  % GetTrajectory get all frames of data for the trial for the specified marker. 
catch % Error check for unlabelled/missing marker
     disp('No T10 trajectory identified - is the marker present and labelled?')
     return
end

idx = T10x<900 & T10x>-900 & T10y>-1500 & T10y<500; % Get index to frames for which T10 is within the selected 3D volume

%Error check for T10 marker not present in capture volume:
if sum(idx)== 0
    disp('No frames within capture volume detected - is the subject within the central test area?')
    return
end
%--------------------
% Knee joint angles
try
[Lcomponents, Le] = GetModelOutput( vicon, Subjects, 'LKneeAngles' );
[Rcomponents, Re] = GetModelOutput( vicon, Subjects, 'RKneeAngles' );
%Returns
            %    components      = numerical(double) NxM matrix where N is the number of components, M is the number of frames
            %    e               = logical array, T/F indication as to whether the data exists for each frame;
catch % Error check for unmodelled data
    disp('No knee angle data found - is the trial modelled?')
    return
end
Le = Le & idx;   Re = Re & idx;  % Only take values that have modelled knee data and are within defined capture volume region           

Lcomponents = [Lcomponents(1,Le); Lcomponents(2,Le) ; Lcomponents(3,Le)]; % Exclude frames with no modelled data for Left and Right
Rcomponents = [Rcomponents(1,Re); Rcomponents(2,Re) ; Rcomponents(3,Re)];

% Get peak frontal and transverse plane knee angles
LPeakIntRot = max(Lcomponents(3,:)) ; LPeakExtRot = min(Lcomponents(3,:)); LRotRange = abs(LPeakIntRot - LPeakExtRot);
LPeakAdd = max(Lcomponents(2,:)) ; LPeakAbd = min(Lcomponents(2,:)); LAbAdRange = abs(LPeakAdd - LPeakAbd);
RPeakIntRot = max(Rcomponents(3,:)) ; RPeakExtRot = min(Rcomponents(3,:)); RRotRange = abs(RPeakIntRot - RPeakExtRot);
RPeakAdd = max(Rcomponents(2,:)) ; RPeakAbd = min(Rcomponents(2,:)); RAbAdRange = abs(RPeakAdd - RPeakAbd);

% ------ Get internal rotation ankle at peak knee flexion and peak knee extension
% Get max and min knee sagittal plane angles and their indexes
[LPkFle, LPkFleIDX] = max(Lcomponents(1,:)); 
[LPkExt, LPkExtIDX] = min(Lcomponents(1,:));
[RPkFle, RPkFleIDX] = max(Rcomponents(1,:));
[RPkExt, RPkExtIDX] = min(Rcomponents(1,:));

% Get knee internal rotation angle at knee max flexion and max extension
LIntRotPkFle = Lcomponents(3,LPkFleIDX);
LIntRotPkExt = Lcomponents(3,LPkExtIDX);
RIntRotPkFle = Rcomponents(3,RPkFleIDX);
RIntRotPkExt = Rcomponents(3,RPkExtIDX);

% Get absolute asymmetry magnitude for internal rotation angle at max fle
% and max ext
 LvRIntRotPkFle = LIntRotPkFle-RIntRotPkFle;
 LvRIntRotPkExt = LIntRotPkExt-RIntRotPkExt;

% Get range of internal rotation from max fle to max ext
LIntRotRange = LIntRotPkFle-LIntRotPkExt ;
RIntRotRange = RIntRotPkFle-RIntRotPkExt;

% figure (1)
subplot(2,1,1)
plot(1:length(Lcomponents), Lcomponents)   
subplot(2,1,2)
plot(1:length(Rcomponents), Rcomponents)   

% disp('--------------------------------------------------------')
% CheckResults = ...
%  [sprintf('Left Knee Peak Rotation Angle (Internal) is %.1f degrees ', LPeakIntRot) ... 
%  sprintf('\nLeft Knee Peak Rotation Angle (External) is %.1f degrees ', LPeakExtRot)...
%  sprintf('\nLeft Knee Rotation Angle Range is %.1f degrees         ', LRotRange)... 
%  sprintf('\n\nRight Knee Peak Rotation Angle (Internal) is %.1f degrees', RPeakIntRot)...
%  sprintf('\nRight Knee Peak Rotation Angle (External) is %.1f degrees', RPeakExtRot)...
%  sprintf('\nRight Knee Rotation Angle Range is %.1f degrees        ', RRotRange)...  
%  sprintf('\n\nLeft Peak Knee Angle (Adduction) is %.1f degrees         ', LPeakAdd)...
%  sprintf('\nLeft Peak Knee Angle (Abduction) is %.1f degrees         ', LPeakAbd)...
%  sprintf('\nLeft Knee Frontal Angle Range is %.1f degrees          ', LAbAdRange)...
%  sprintf('\n\nRight Peak Knee Angle (Adduction) is %.1f degrees        ',RPeakAdd)...
%  sprintf('\nRight Peak Knee Angle (Abduction) is %.1f degrees        ', RPeakAbd)...
%  sprintf('\nRight Knee Frontal Angle Range is %.1f degrees         ', RAbAdRange)...
%  ]
% disp('--------------------------------------------------------')

 
 disp('--------------------------------------------------------')
CheckResults2 = ...
 [sprintf('Knee Rotation Asymmetry at max flexion is %.1f degrees (positive value means Left more internally-rotated than Right)  ', LvRIntRotPkFle) ... 
 sprintf('\nKnee Rotation Asymmetry at max extension is %.1f degrees (positive value means Left more internally-rotated than Right) ', LvRIntRotPkExt)...
 sprintf('\nBOTH MUST BE BETWEEN %.1f AND %.1f DEGREES', -10, 10)...
 sprintf('\n\nLeft Knee Rotation Angle Range (flex-ext) is %.1f degrees', LIntRotRange)...
 sprintf('\nRight Knee Rotation Angle Range (flex-ext) is %.1f degrees', RIntRotRange)...
 sprintf('\nBOTH MUST BE LESS THAN %.1f DEGREES', 20)...
 ]
disp('--------------------------------------------------------')
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
            
            % %%--------------------
% % Write to ENF
% 
% [ path, name ] = GetTrialName( vicon );
%             % GetTrialName retrieves the name and path for the loaded trial            
%             %
%             % Returns
%             %    path = string, path to the trial on disk
%             %    name = string, name of the trial
% 
% path = path(1:end-1); %remove trailing backslash
% files = searchFolder4Files(path);
% ispresent = cellfun(@(s) ~isempty(strfind(s, name)), files) & cellfun(@(s) ~isempty(strfind(s, '.enf')), files); % get index to enf file of current trial in files list
% enfTrial = files(ispresent); % name of trial ENF file
% enfTrial = char(enfTrial);
% enftool = EnfToolbox;
% 
% enffile = enfTrial;
% term = 'hi loser';
% 
% 
%             
%           % open enffile
%             fileID = fopen(enffile,'r');
%           % read file line by line
%             criteria = true; rep = 1; newENF = '';
%             while criteria
%                 newENF{rep,1} = fgets(fileID);
%                 if newENF{rep} == -1
%                     criteria = false;
%                     newENF = newENF(1:rep-1,1);
%                 end
%                 if rep == 1
%                     seperator = newENF{rep,1}(end-1:end);
%                 end
%                     rep = rep + 1;
%             end
%           % close enffile                
%             fclose(fileID);
%           % check if enf has already a note
%             if sum(cell2mat(strfind(newENF,'NOTES='))) == 1
%                 count = 1;
%                 while count <= size(newENF,1)
%                     if ~isempty(strfind(newENF{count,1},'NOTES='))
%                         newENF{count,1} = ['NOTES=',term,seperator];
%                     end
%                     count = count+1;
%                 end
%             elseif sum(cell2mat(strfind(newENF,'NOTES='))) == 0  
%                 newENF = enftool.addNoteToENFfile(enffile,newENF,term,seperator);
%             elseif sum(cell2mat(strfind(newENF,'NOTES='))) >= 2
%                 count = 1; idx = zeros(size(newENF));
%                 while count <= size(newENF,1)
%                     if ~isempty(strfind(newENF{count,1},'NOTES='))
%                         idx(count) = 1;
%                     end
%                     count = count+1;
%                 end
%                 newENF = newENF(idx==0);
%                 newENF = enftool.addNoteToENFfile(enffile,newENF,term,seperator);
%             end
%             
% 
%             
%           % overwrite enffile with new info                
%             txt = fopen(enffile,'w');   
%             for n = 1:size(newENF,1)
%                 fprintf(txt,'%s',newENF{n,:});
%             end
%             fclose(txt);

