%This code is designed to analysis selected repeated jump trials from vicon
%3D files (.Enf), once files are selected the code will open each trial on
%vicon, calculate contact time, jump height and reactive strength index of
%each each jump from trial, then write the results into an excel file.
%Authors - Yunus Kaghembe and Leigh James Ryan
%Created 23/03/2021

clear
clc

try
    vicon= ViconNexus;
catch
    error = errordlg('Please ensure that Vicon is open before executing code', 'Unable to run code');
    uiwait(error)
    return
end
msg=msgbox("Select RJT or DJ trial with injured leg last to run code properly");    
uiwait(msg)

[Fname, Pname] = uigetfile('.Enf', 'Select RJT and DJ Files', 'MultiSelect', 'on');

CF = contains(Fname, 'RJT');
CF = nonzeros(CF);
CF=CF';
              
DF = contains(Fname, 'DJ');
DF= nonzeros(DF);
DF=DF';

       
for p = 1:length(CF)
    if CF(1,p) < 1
        filecheck = msgbox('Please check Fname and ensure that correct RJT files are selected');
        uiwait(filecheck)
        return
    end
end
for p = 1:length(DF)
    if DF(1,p) < 1
        filecheck = msgbox('Please check Fname and ensure that correct DJ files are selected');
        uiwait(filecheck)
        return
    end
end

if ischar(Fname)
    y = 1;
elseif ~ischar(Fname)
    y = size(Fname, 2);
end
% Select injuried leg
[injside,~] = listdlg('ListString', {'Left', 'Right'},'SelectionMode','single', 'PromptString', 'Select the injured side',...
         'ListSize', [150,50], 'CancelString', 'Quit'); 
     
% Initializing DJ variables
RSI3=[];
JumpHeightDJ3=[];
ContactTime3=[];
RSISLDJR=[];
JumpHeightSLDJR=[];
ContactTimeSLDJR=[];
RSISLDJL=[];
JumpHeightSLDJL=[];
ContactTimeSLDJL=[];

% Initializing Joint Moment and Stiffness Variables
DLDJ_ARmom=[];
DLDJ_ALmom=[];
DLDJ_KRmom=[];
DLDJ_KLmom=[];
DLDJ_HRmom=[];
DLDJ_HLmom=[];
SLDJR_ARmom=[];
SLDJR_KRmom=[];
SLDJR_HRmom=[];
SLDJL_ALmom=[];
SLDJL_KLmom=[];
SLDJL_HLmom=[];

AnkleRx_SLRJT=[];
KneeRx_SLRJT=[];
HipRx_SLRJT=[];
AnkleLx_SLRJT=[];
KneeLx_SLRJT=[];
HipLx_SLRJT=[];
AnkleRx_DLRJT=[];
KneeRx_DLRJT=[];
HipRx_DLRJT=[];
AnkleLx_DLRJT=[];
KneeLx_DLRJT=[];
HipLx_DLRJT=[];

DLDJ_ARstiff=[];
DLDJ_ALstiff=[];
DLDJ_KRstiff=[];
DLDJ_KLstiff=[]; 
DLDJ_HRstiff=[]; 
DLDJ_HLstiff=[];
SLDJR_ARstiff=[];
SLDJR_KRstiff=[];
SLDJR_HRstiff=[];
SLDJL_ALstiff=[];
SLDJL_KLstiff=[];
SLDJL_HLstiff=[];
COM_Stiff_DJ=[];
COM_Stiff_SLDJR=[];
COM_Stiff_SLDJL=[];
DJ_Data_L=[];
DJ_Data_R=[];
RJT_Data_L=[];
RJT_Data_R=[];
RAnkleJointStiff_SLRJT=[];
RKneeJointStiff_SLRJT=[];
RHipJointStiff_SLRJT=[];
LAnkleJointStiff_SLRJT=[];
LKneeJointStiff_SLRJT=[];
LHipJointStiff_SLRJT=[];

RAnkleJointStiffness_DLRJT=[];
RKneeJointStiffness_DLRJT=[];
RHipJointStiffness_DLRJT=[];
LAnkleJointStiffness_DLRJT=[];
LKneeJointStiffness_DLRJT=[];
LHipJointStiffness_DLRJT=[];

for n = 1:y
    %% Opening the right trial and getting the required information i.e. Bodymass etc...
    if ischar(Fname)
        answer = questdlg('Only one trial has been selected', 'Trials Selected', 'Continue' ,'Cancel', 'Cancel');
        if strcmp(answer,'Cancel')
            return
        elseif strcmp(answer,'Continue')
            idx = find(Fname == '.', 1, 'first');
            Tname = Fname(1:idx-1);
        end
    elseif ~ischar(Fname)
        x = Fname{1,n};
        idx = find(x == '.', 1, 'first');
        Tname = x(1:idx-1); %Name of Trial
    end
    
    Trialpath = strcat(Pname, Tname); %Trial path
    vicon.OpenTrial(Trialpath, 20); % opening trial
    
    %getting frames
    [sFrame,eFrame] = vicon.GetTrialRegionOfInterest;
    Fs = vicon.GetFrameRate; %Sampling Frequency
    
    %Get Subjects name
    SubjectName = vicon.GetSubjectNames;
    SubjectName = char (SubjectName(1,1)); %converting the subject name into a character
    
    %Getting Subject Parameters
    SubjectParam = vicon.GetSubjectParamNames(SubjectName);
    BM = vicon.GetSubjectParam(SubjectName,SubjectParam{1,1}); %Bodymass
    BW = BM*9.81; %Bodyweight
    %Getting Model Parameters
    COMxyz=vicon.GetModelOutput(SubjectName,'CentreOfMass');
    RKneeAngle=vicon.GetModelOutput(SubjectName,'RKneeAngles');
    LKneeAngle=vicon.GetModelOutput(SubjectName,'LKneeAngles');
    RAnkleAngle=vicon.GetModelOutput(SubjectName,'RAnkleAngles');
    LAnkleAngle=vicon.GetModelOutput(SubjectName,'LAnkleAngles')
    RHipAngle=vicon.GetModelOutput(SubjectName,'RHipAngles');
    LHipAngle=vicon.GetModelOutput(SubjectName,'LHipAngles');
    KneeR= vicon.GetModelOutput(SubjectName,'RKneeMoment');
    KneeL= vicon.GetModelOutput(SubjectName,'LKneeMoment');
    AnkleR= vicon.GetModelOutput(SubjectName,'RAnkleMoment');
    AnkleL= vicon.GetModelOutput(SubjectName,'LAnkleMoment');
    HipR= vicon.GetModelOutput(SubjectName,'RHipMoment');
    HipL= vicon.GetModelOutput(SubjectName,'LHipMoment');
    %% importing analog force data
    
    %Device ID
    ID = vicon.GetDeviceIDs;
    DeviceID = cell(length(ID),1);
    for x = 1:length(DeviceID)
        DeviceID{x} = ID(x);
    end
    %Channels
    Channels = {'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};
    
    %Preallocating
    %DeviceID = cell(length(DeviceName),1);
    DeviceType = cell(length(DeviceID),1);
    DeviceRate = cell(length(DeviceID),1);
    DeviceOutputID = cell(length(DeviceID),1);
    ChannelID = cell(length(Channels),1);
    
    for i = 1:length(DeviceID)
        %DeviceID{i} = vicon.GetDeviceIDFromName(DeviceName{i});
        [~,DeviceType{i}, DeviceRate{i}, DeviceOutputID{i}] = vicon.GetDeviceDetails(DeviceID{i});
    end
    
    for i = 1:length(DeviceID)
        if strcmp(DeviceType{i},'ForcePlate');
            for j = 1:length(Channels)
                ChannelID{j} = vicon.GetDeviceChannelIDFromName(DeviceID{i}, ceil(double(j)/3), Channels{j});
                % 2nd term is Device OutputID (1 = force, 2 = moment, 3 = CoP)
                ForcePlateRaw.DeviceID{i}.Channels{j} = vicon.GetDeviceChannel(DeviceID{i}, ceil(double(j)/3), ChannelID{j});
                %ForcePlateCropped.DeviceID{i}.Channels{j} = ForcePlateRaw.DeviceID{i}.Channels{j}(sFrame*(DeviceRate{i}/Fs):eFrame*(DeviceRate{i}/Fs));
            end
        end
    end
    
    Fz2 = ForcePlateRaw.DeviceID{1,1}.Channels{1,3} * -1;
    Fz1 = ForcePlateRaw.DeviceID{1,2}.Channels{1,3} * -1;
    
    VerticalForce = Fz1 + Fz2; %Combine both forceplates as one forceplate
    VerticalForce = VerticalForce';
    %TrialForce{n} = VerticalForce;
    
    %donwsample forcedata into 200 Hz similar to kinematic data
    VerticalForce = decimate(VerticalForce, 5);
    VerticalForce = VerticalForce(sFrame:eFrame);
    
    %% Classifying Events
    [FS, FO]= getevents(VerticalForce, BW);
    k = 1:length(VerticalForce);
    plot(VerticalForce)
    title(Tname)
    hold on
    plot(k(FS), VerticalForce(FS),'x')
    plot(k(FO), VerticalForce(FO),'o')
    hold off
    
    pause;
    close
    %% Calculating RJT and DJ Variables and Selecting top five RJT jumps
    
    [ContactTime,JumpHeightDJ,JumpHeightv1,RSI,COMv0,COM_Stiff,COM_Stiff_RJT,COMv0_RJT] = RJTanalysis(FS,FO,COMxyz,VerticalForce,BW,BM);
    
   
    
    
    
   %Conjugating DJ jump variables
    
    if length(FO)<2
        r=RSI;
        dj=JumpHeightDJ;
        c=ContactTime;
        if contains(Tname,'SLDJRight')==1
            RSISLDJR = [RSISLDJR,r]
            JumpHeightSLDJR = [JumpHeightSLDJR,dj]
            ContactTimeSLDJR= [ContactTimeSLDJR,c]
        elseif contains(Tname,'SLDJLeft')==1
            RSISLDJL = [RSISLDJL,r]
            JumpHeightSLDJL = [JumpHeightSLDJL,dj]
            ContactTimeSLDJL = [ContactTimeSLDJL,c]
        else       
            RSI3=[RSI3,r];
            JumpHeightDJ3=[JumpHeightDJ3,dj];
            ContactTime3=[ContactTime3,c];
        end
    end
        
        
    if length(FO)>2
    % top five RJT jumps
    [jumps, index] = sort(RSI); %rearrange the jumps in order of lowest to highest RSI
    newindex = index(:,end-4:end);
    
       for v = 1:5
          jh = JumpHeightv1(:,newindex(v));
          JH5(:,v) = jh;
        
          ct = ContactTime(:, newindex(v));
          CT5(:,v) = ct;
         
          rsi = RSI(:, newindex(v));
          RSI5(:,v) = rsi;
       end
    
    
    %mean of the 5 best RJT jumps
    
      AveCT = round(mean(CT5),3);
      AveJH = round(mean(JH5),3);
      AveRSI = round(mean(RSI5),3);
      
    end
    
    %% Writing Results into Excel file
   
    %header1 ={ 'DJ and RJT jump variables'};
    header2 = {'JumpHeight (m)' 'Contact Time (s)' 'RSI'}';
    header3 = {'Best 5'};
    
    if length(FO)>2
    for q = 1:length(JumpHeightv1);
        NO = num2str(q);
        hd = ['Jump ' NO];
        header4(:,q) = {hd};
    end
  
  
    
    for a =1:5
       hd2 = header4(:,newindex(:,a));
    header5(:,a) = hd2; 
    end 
   
   header6 = {'Average'};
    
    Data =[JumpHeightv1; ContactTime; RSI];
    Data2 = [JH5; CT5; RSI5];
    Data3 = [AveJH; AveCT; AveRSI];
   
   header7 = {'Injuried'};
    
    filename = strcat(SubjectName,'RJT');
     if injside == 1 
        xlswrite(filename,header7,"SLRJTLeft","A1")
     else 
        xlswrite(filename,header7,"SLRJTRight","A1")
     end
   header8= {'LSI Values'}; 
          
         
   
    xlswrite(filename,header4,Tname,'B1')
    xlswrite(filename,header2,Tname,'A2')
    xlswrite(filename,header3,Tname,'A7')
    xlswrite(filename,header5,Tname,'B7')
    %xlswrite(filename,header1,Tname,'A8')
    xlswrite(filename,Data,Tname,'B2')
    xlswrite(filename,Data2,Tname,'B8')
    xlswrite(filename,header6,Tname,'H7')
    xlswrite(filename,Data3,Tname,'H8')
    
    if injside==1
    xlswrite(filename,header7,Tname,'I7')
    else
    xlswrite(filename,header7,Tname,'I7')
    end
   %Calculating LSI
    if contains(Tname,'SLRJTLeft')==1
        DataLeft=Data3;
        RJT_Data_L=[RJT_Data_L,DataLeft];
    elseif contains(Tname,'SLRJTRight')==1
        DataRight=Data3;
        RJT_Data_R=[RJT_Data_R,DataRight];
    end
    
    
    
    
    
    
      
    %if contains(Tname,"Left")== 1
     %cells2 = xlsread(filename,Tname,"H8:H10");
      
    %elseif contains(Tname,"Right")== 1 
     %cells3 = xlsread(filename,Tname,"H8:H10");
    
    %else
       %cells2 = xlsread(filename,Tname,"H8:H10");
       %cells3 = xlsread(filename,Tname,"H8:H10");
     
   
    %end 
    
    %if injside==1 && contains(Tname,"Left")==1
      %LSI=(cells2./cells3)*100;
      %LSI1= round(LSI,2);
      %xlswrite(filename,LSI1,Tname,'I8')
    %elseif  injside==2 && contains(Tname,"Right")==1 
      %LSI=(cells3./cells2)*100;
      %LSI1= round(LSI,2);
      %xlswrite(filename,LSI1,Tname,'I8')
    
    
    %end
    
    
    
  elseif length(FO)<2
       header1 ={ 'DJ jump variables'};
       header2 = {'JumpHeight (m)' 'Contact Time (s)' 'RSI'}';
       header6 = {'Average'};
        if contains(Tname,'DLDJ')==1 && numel(RSI3)>2
             filename = strcat(SubjectName,'RJT')
             Data =[JumpHeightDJ3; ContactTime3; RSI3];
             JumpHeightDJ3_aver=mean(JumpHeightDJ3);
             ContactTime3_aver=mean(ContactTime3);
             RSI3_aver=mean(RSI3);
             AVER=[JumpHeightDJ3_aver;ContactTime3_aver;RSI3_aver]
             xlswrite(filename,header2,Tname,"A1:A3")
             xlswrite(filename,Data,Tname,"B1:D3");
             xlswrite(filename,header6,Tname,"C4");
             xlswrite(filename,AVER,Tname,"C5:C7");
        
        elseif contains(Tname,'SLDJRight')==1 && numel(RSISLDJR)>2
             header6 = {'Average'};
             filename = strcat(SubjectName,'RJT')
             Data =[JumpHeightSLDJR; ContactTimeSLDJR; RSISLDJR];
             JumpHeightSLDJR_aver=mean(JumpHeightSLDJR);
             ContactTimeSLDJR_aver=mean(ContactTimeSLDJR);
             RSISLDJR_aver=mean(RSISLDJR);
             AVER_R=[JumpHeightSLDJR_aver;ContactTimeSLDJR_aver;RSISLDJR_aver]
             xlswrite(filename,header2,Tname,"A1:A3")
             xlswrite(filename,Data,Tname,"B1:D3");
             xlswrite(filename,header6,Tname,"C4");
             xlswrite(filename,AVER_R,Tname,"C5:C7");
        
        elseif contains(Tname, 'SLDJLeft')==1 && numel(RSISLDJL)>2
             header6 = {'Average'};
             filename = strcat(SubjectName,'RJT')
             Data =[JumpHeightSLDJL; ContactTimeSLDJL; RSISLDJL];
             JumpHeightSLDJL_aver=mean(JumpHeightSLDJL);
             ContactTimeSLDJL_aver=mean(ContactTimeSLDJL);
             RSISLDJL_aver=mean(RSISLDJL);
             AVER_L=[JumpHeightSLDJL_aver;ContactTimeSLDJL_aver;RSISLDJL_aver]
             xlswrite(filename,header2,Tname,"A1:A3")
             xlswrite(filename,Data,Tname,"B1:D3");
             xlswrite(filename,header6,Tname,"C4");
             xlswrite(filename,AVER_L,Tname,"C5:C7");
       end      
   else
       continue
   end    
RKneeAnglexStance=[];
LKneeAnglexStance=[];
RAnkleAnglexStance=[];
LAnkleAnglexStance=[];
RHipAnglexStance=[];
LHipAnglexStance=[];

AnkleRxStance=[];
AnkleLxStance=[];
KneeRxStance=[];
KneeLxStance=[];
HipRxStance=[];
HipLxStance=[];
 %% Plotting time-normalised Sagittal-plane joint moment
 %Getting sagittal plane angles and moments
 RKneeAnglex=RKneeAngle(1,:);
 LKneeAnglex=LKneeAngle(1,:);
 RAnkleAnglex=RAnkleAngle(1,:);
 LAnkleAnglex=LAnkleAngle(1,:);
 RHipAnglex= RHipAngle(1,:);
 LHipAnglex=LHipAngle(1,:);
 
 AnkleRx=AnkleR(1,:);
 AnkleLx=AnkleL(1,:);
 KneeRx=KneeR(1,:);
 KneeLx =KneeL(1,:);
 HipRx =HipR(1,:);
 HipLx = HipL(1,:);
 
 %finding Stance phase
  RKneeAnglexS=RKneeAnglex(FS(1):FO(1));
  LKneeAnglexS=LKneeAnglex(FS(1):FO(1));
  RAnkleAnglexS=RAnkleAnglex(FS(1):FO(1));
  LAnkleAnglexS=LAnkleAnglex(FS(1):FO(1));
  RHipAnglexS= RHipAnglex(FS(1):FO(1));
  LHipAnglexS=LHipAnglex(FS(1):FO(1));
    
  AnkleRxS=AnkleRx(FS(1):FO(1));
  AnkleLxS=AnkleLx(FS(1):FO(1));
  KneeRxS=KneeRx(FS(1):FO(1));
  KneeLxS=KneeLx(FS(1):FO(1));
  HipRxS=HipRx(FS(1):FO(1));
  HipLxS=HipLx(FS(1):FO(1));
  
 if contains(Tname,'SLRJT') || contains(Tname,'RJT')
 [num,txt]=xlsread(filename1 ,Tname,'B7:F7');
 txt=erase(txt,'Jump');
 txt=str2double(txt);
  for n=1:length(txt)
    RKneeAnglexS=RKneeAnglex(FS(txt(n)):FO(txt(n)+1));
    LKneeAnglexS=LKneeAnglex(FS(txt(n)):FO(txt(n)+1));
    RAnkleAnglexS=RAnkleAnglex(FS(txt(n)):FO(txt(n)+1));
    LAnkleAnglexS=LAnkleAnglex(FS(txt(n)):FO(txt(n)+1));
    RHipAnglexS= RHipAnglex(FS(txt(n)):FO(txt(n)+1))*1/10;
    LHipAnglexS=LHipAnglex(FS(txt(n)):FO(txt(n)+1));
    
    
  
    AnkleRxS=AnkleRx(FS(txt(n)):FO(txt(n)+1))*1/1000;
    AnkleLxS=AnkleLx(FS(txt(n)):FO(txt(n)+1))*1/1000;
    KneeRxS=KneeRx(FS(txt(n)):FO(txt(n)+1))*1/1000;
    KneeLxS=KneeLx(FS(txt(n)):FO(txt(n)+1))*1/1000;
    HipRxS=HipRx(FS(txt(n)):FO(txt(n)+1))*1/1000;
    HipLxS=HipLx(FS(txt(n)):FO(txt(n)+1))*1/1000;
    
    
  end
  RKneeAnglexStance=[RKneeAnglexStance;RKneeAnglexS];
  LKneeAnglexStance=[LKneeAnglexStance;LKneeAnglexS];
  RAnkleAnglexStance=[RAnkleAnglexStance;RAnkleAnglexS];
  LAnkleAnglexStance=[LAnkleAnglexStance;LAnkleAnglexS];
  RHipAnglexStance=[RHipAnglexStance;RHipAnglexS];
  LHipAnglexStance=[LHipAnglexStance;LHipAnglexS];
  
  AnkleRxStance=[AnkleRxStance;AnkleRxS];
  AnkleLxStance=[AnkleLxStance;AnkleLxS];
  KneeRxStance=[KneeRxStance;KneeRxS];
  KneeLxStance=[KneeLxStance;KneeLxS];
  HipRxStance=[HipRxStance;HipRxS];
  HipLxStance=[HipLxStance;HipLxS];
  
  
  Aver_RKneeAnglexStance=mean(RKneeAnglexStance);
  Aver_LKneeAnglexStance=mean(LKneeAnglexStance);
  Aver_RAnkleAnglexStance=mean(RAnkleAnglexStance);
  Aver_LAnkleAnglexStance=mean(LAnkleAnglexStance);
  Aver_RHipAnglexStance=mean(RHipAnglexStance);
  Aver_LHipAnglexStance=mean(LHipAnglexStance);
  
  
  
  Aver_AnkleRxStance=mean(AnkleRxStance);
  Aver_AnkleLxStance=mean(AnkleLxStance);
  Aver_KneeRxStance=mean(KneeRxStance);
  Aver_KneeLxStance=mean(KneeLxStance);
  Aver_HipRxStance=mean(HipRxStance);
  Aver_HipLxStance=mean(HipLxStance);
 end
 
 
 
 
 
 heading= {'LAnkle Joint Moment'};
 heading1= {'LKnee Joint Moment'};
 heading2= {'LHip Joint Moment'};
 heading3= {'RAnkle Joint Moment'};
 heading4= {'RKnee Joint Moment'};
 heading5= {'RHip Joint Moment'};
 
 heading6={'Aver LAnkle Joint Moment'};
 heading7={'Aver LKnee Joint Moment'};
 heading8={'Aver LHip Joint Moment'};
 
 heading9={'Aver RAnkle Joint Moment'};
 heading10={'Aver RKnee Joint Moment'};
 heading11={'Aver RHip Joint Moment'};
 
 heading12= {'LAnkle Joint Stiffness'};
 heading13= {'LKnee Joint Stiffness'};
 heading14= {'LHip Joint Stiffness'};
 heading15= {'RAnkle Joint Stiffness'};
 heading16= {'RKnee Joint Stiffness'};
 heading17= {'RHip Joint Stiffness'};
 heading18={'Average'}
 heading19= {'COM Stiffness'};
 
 
 
 filename = strcat(SubjectName,'RJT');
  if contains(Tname,'DLDJ')==1
      AnkleRxDJ=interp1(1:numel(AnkleRxS),AnkleRxS, linspace(1, numel(AnkleRxS), 101))*1/1000;
      AnkleLxDJ=interp1(1:numel(AnkleLxS),AnkleLxS, linspace(1, numel(AnkleLxS), 101))*1/1000;
      KneeRxDJ=interp1(1:numel(KneeRxS),KneeRxS, linspace(1, numel(KneeRxS), 101))*1/1000;
      KneeLxDJ=interp1(1:numel(KneeLxS),KneeLxS, linspace(1, numel(KneeLxS), 101))*1/1000;
      HipRxDJ=interp1(1:numel(HipRxS),HipRxS, linspace(1, numel(HipRxS), 101))*1/1000;
      HipLxDJ=interp1(1:numel(HipLxS),HipLxS, linspace(1, numel(HipLxS), 101))*1/1000;
      
      DLDJ_ARmom=[DLDJ_ARmom;AnkleRxDJ];
      DLDJ_ALmom=[DLDJ_ALmom; AnkleLxDJ];
      DLDJ_KRmom=[DLDJ_KRmom; KneeRxDJ];
      DLDJ_KLmom=[DLDJ_KLmom; KneeLxDJ]; 
      DLDJ_HRmom=[DLDJ_HRmom; HipRxDJ]; 
      DLDJ_HLmom=[DLDJ_HLmom; HipLxDJ];
      
%********************DLDJ Joint Stiffness****************************
        %Ankle Joint Stiffness
        %Right
       RAnkleAnglexIC=RAnkleAnglex(FS);
       RAnkleAnglexPeakFlex = RAnkleAnglexS(COMv0);
       Delta_RAnkleAngle=RAnkleAnglexPeakFlex-RAnkleAnglexIC
       Delta_AnkleRxS=AnkleRxS(COMv0)*1/1000;
       RAnkleJointStiff=Delta_AnkleRxS/Delta_RAnkleAngle;
        %Left
       LAnkleAnglexIC=LAnkleAnglex(FS);
       LAnkleAnglexPeakFlex = LAnkleAnglexS(COMv0);
       Delta_LAnkleAngle=LAnkleAnglexPeakFlex-LAnkleAnglexIC
       Delta_AnkleLxS=AnkleLxS(COMv0)*1/1000;
       LAnkleJointStiff=Delta_AnkleLxS/Delta_LAnkleAngle;
            
        %Knee Joint Stiffness
        %Right
       RKneeAnglexIC=RKneeAnglex(FS);
       RKneeAnglexPeakFlex = RKneeAnglexS(COMv0);
       Delta_RKneeAngle=RKneeAnglexPeakFlex-RKneeAnglexIC
       Delta_KneeRxS=KneeRxS(COMv0)*1/1000;
       RKneeJointStiff=Delta_KneeRxS/Delta_RKneeAngle;
        %Left
       LKneeAnglexIC=LKneeAnglex(FS);
       LKneeAnglexPeakFlex = LKneeAnglexS(COMv0);
       Delta_LKneeAngle=LKneeAnglexPeakFlex-LKneeAnglexIC;
       Delta_KneeLxS=KneeLxS(COMv0)*1/1000;
       LKneeJointStiff=Delta_KneeLxS/Delta_LKneeAngle;
       
        %Hip Joint Stiffness
        %Right
       RHipAnglexIC=RHipAnglex(FS);
       RHipAnglexPeakFlex=RHipAnglexS(COMv0);
       Delta_RHipAngle=RHipAnglexPeakFlex-RHipAnglexIC;
       Delta_HipRxS=HipRxS(COMv0)*1/1000;
       RHipJointStiff=Delta_HipRxS/Delta_RHipAngle;
        %Left
       LHipAnglexIC=LHipAnglex(FS);
       LHipAnglexPeakFlex = LHipAnglexS(COMv0);
       Delta_LHipAngle=LHipAnglexPeakFlex-LHipAnglexIC;
       Delta_HipLxS=HipLxS(COMv0)*1/1000;
       LHipJointStiff=Delta_HipLxS/Delta_LHipAngle; 
        
      DLDJ_ARstiff=[DLDJ_ARstiff,RAnkleJointStiff];
      DLDJ_ALstiff=[DLDJ_ALstiff,LAnkleJointStiff];
      DLDJ_KRstiff=[DLDJ_KRstiff,RKneeJointStiff];
      DLDJ_KLstiff=[DLDJ_KLstiff,LKneeJointStiff]; 
      DLDJ_HRstiff=[DLDJ_HRstiff,RHipJointStiff]; 
      DLDJ_HLstiff=[DLDJ_HLstiff,LHipJointStiff];
      %Calculating COM Stiffness
      COM_Stiff_DJ=[COM_Stiff_DJ,COM_Stiff];
      
      
      if contains(Tname,'DLDJ3')==1
         DLDJ_ALmom_aver= mean(DLDJ_ALmom); 
         DLDJ_ARmom_aver= mean(DLDJ_ARmom); 
         DLDJ_KLmom_aver= mean(DLDJ_KLmom);
         DLDJ_KRmom_aver= mean(DLDJ_KRmom);
         DLDJ_HLmom_aver= mean(DLDJ_HLmom);
         DLDJ_HRmom_aver= mean(DLDJ_HRmom);
         
         DLDJ_ARstiff_aver= mean(DLDJ_ARstiff);
         DLDJ_ALstiff_aver= mean(DLDJ_ALstiff);
         DLDJ_KRstiff_aver= mean(DLDJ_KRstiff);
         DLDJ_KLstiff_aver= mean(DLDJ_KLstiff);
         DLDJ_HRstiff_aver= mean(DLDJ_HRstiff);
         DLDJ_HLstiff_aver= mean(DLDJ_HLstiff);
         
         xlswrite(filename,heading,Tname,'A12');
         xlswrite(filename,heading6,Tname, 'A15');
         xlswrite(filename,heading1,Tname, 'A16');
         xlswrite(filename,heading7,Tname, 'A19');
         xlswrite(filename,heading2,Tname, 'A20');
         xlswrite(filename,heading8,Tname, 'A23');
         xlswrite(filename,heading3,Tname, 'A24');
         xlswrite(filename,heading9,Tname, 'A27');
         xlswrite(filename,heading4,Tname, 'A28');
         xlswrite(filename,heading10,Tname, 'A31');
         xlswrite(filename,heading5,Tname, 'A32');
         xlswrite(filename,heading11,Tname, 'A35');
         
         xlswrite(filename,heading12,Tname,'F2');
         xlswrite(filename,heading13,Tname, 'F3');
         xlswrite(filename,heading14,Tname, 'F4');
         xlswrite(filename,heading15,Tname, 'F5');
         xlswrite(filename,heading16,Tname, 'F6');
         xlswrite(filename,heading17,Tname, 'F7');
         xlswrite(filename,heading18,Tname, 'J1');
         
         xlswrite(filename,DLDJ_ALmom,Tname, 'B12:CX14');
         xlswrite(filename,DLDJ_ALmom_aver,Tname, 'B15:CX15');
         xlswrite(filename,DLDJ_KLmom,Tname, 'B16:CX18');
         xlswrite(filename,DLDJ_KLmom_aver,Tname, 'B19:CX19');
         xlswrite(filename,DLDJ_HLmom,Tname, 'B20:CX22');
         xlswrite(filename,DLDJ_HLmom_aver,Tname, 'B23:CX23');
         xlswrite(filename,DLDJ_ARmom,Tname, 'B24:CX26');
         xlswrite(filename,DLDJ_ARmom_aver,Tname, 'B27:CX27');
         xlswrite(filename,DLDJ_KRmom,Tname, 'B28:CX30');
         xlswrite(filename,DLDJ_KRmom_aver,Tname, 'B31:CX31');
         xlswrite(filename,DLDJ_HRmom,Tname, 'B32:CX34');
         xlswrite(filename,DLDJ_HRmom_aver,Tname, 'B35:CX35');
         
         xlswrite(filename,DLDJ_ALstiff,Tname, 'G2:I2');
         xlswrite(filename,DLDJ_ALstiff_aver,Tname, 'J2');
         xlswrite(filename,DLDJ_KLstiff,Tname, 'G3');
         xlswrite(filename,DLDJ_KLstiff_aver,Tname, 'J3');
         xlswrite(filename,DLDJ_HLstiff,Tname, 'G4');
         xlswrite(filename,DLDJ_HLstiff_aver,Tname, 'J4');
         xlswrite(filename,DLDJ_ARstiff,Tname, 'G5');
         xlswrite(filename,DLDJ_ARstiff_aver,Tname, 'J5');
         xlswrite(filename,DLDJ_KRstiff,Tname, 'G6');
         xlswrite(filename,DLDJ_KRstiff_aver,Tname, 'J6');
         xlswrite(filename,DLDJ_HRstiff,Tname, 'G7');
         xlswrite(filename,DLDJ_HRstiff_aver,Tname, 'J7');
         
         %COM Stiffness
         xlswrite(filename,heading19,Tname, 'L1');
         xlswrite(filename,COM_Stiff_DJ,Tname, 'L2');
         xlswrite(filename,heading18,Tname, 'O1');
         COM_Stiff_DJ_aver=mean(COM_Stiff_DJ);
         xlswrite(filename,COM_Stiff_DJ_aver,Tname, 'O2');
      end
      
  elseif contains(Tname, 'SLDJRight')==1
      AnkleRxSLDJR=interp1(1:numel(AnkleRxS),AnkleRxS, linspace(1, numel(AnkleRxS), 101))*1/1000;
      KneeRxSLDJR=interp1(1:numel(KneeRxS),KneeRxS, linspace(1, numel(KneeRxS), 101))*1/1000;    
      HipRxSLDJR=interp1(1:numel(HipRxS),HipRxS, linspace(1, numel(HipRxS), 101))*1/1000;
      
      SLDJR_ARmom=[SLDJR_ARmom; AnkleRxSLDJR];
      
      SLDJR_KRmom=[SLDJR_KRmom; KneeRxSLDJR];
     
      SLDJR_HRmom=[SLDJR_HRmom; HipRxSLDJR];
      
      
      
%*****************SLDJRight Joint Stiffness***************       
        %Ankle Joint Stiffness
        %Right
       RAnkleAnglexIC=RAnkleAnglex(FS);
       RAnkleAnglexPeakFlex =RAnkleAnglexS(COMv0);
       Delta_RAnkleAngle=RAnkleAnglexPeakFlex-RAnkleAnglexIC;
       Delta_AnkleRxS=AnkleRxS(COMv0)*1/1000;
       RAnkleJointStiff=Delta_AnkleRxS/Delta_RAnkleAngle;
        
       
        %Knee Joint Stiffness
        %Right
       RKneeAnglexIC=RKneeAnglex(FS);
       RKneeAnglexPeakFlex = RKneeAnglexS(COMv0);
       Delta_RKneeAngle=RKneeAnglexPeakFlex-RKneeAnglexIC;
       Delta_KneeRxS=KneeRxS(COMv0)*1/1000;
       RKneeJointStiff=Delta_KneeRxS/Delta_RKneeAngle;
        
       
        %Hip Joint Stiffness
        %Right
       RHipAnglexIC=RHipAnglex(FS);
       RHipAnglexPeakFlex =RHipAnglexS(COMv0);
       Delta_RHipAngle=RHipAnglexPeakFlex-RHipAnglexIC;
       Delta_HipRxS=HipRxS(COMv0)*1/1000;
       RHipJointStiff=Delta_HipRxS/Delta_RHipAngle;
       
       
       SLDJR_ARstiff=[SLDJR_ARstiff,RAnkleJointStiff];
      
       SLDJR_KRstiff=[SLDJR_KRstiff,RKneeJointStiff];
     
       SLDJR_HRstiff=[SLDJR_HRstiff,RHipJointStiff];
       
       %Calculating COM Stiffness
       COM_Stiff_SLDJR=[COM_Stiff_SLDJR,COM_Stiff];
       
       
      
      if contains(Tname,'SLDJRight3')==1
          SLDJR_ARmom_aver=mean(SLDJR_ARmom);
          SLDJR_KRmom_aver=mean(SLDJR_KRmom);
          SLDJR_HRmom_aver=mean(SLDJR_HRmom);
          
          SLDJR_ARstiff_aver= mean(SLDJR_ARstiff);        
          SLDJR_KRstiff_aver= mean(SLDJR_KRstiff);
          SLDJR_HRstiff_aver= mean(SLDJR_HRstiff);
         
         
          xlswrite(filename,heading3,Tname, 'A24');
          xlswrite(filename,heading9,Tname, 'A27');
          xlswrite(filename,heading4,Tname, 'A28');
          xlswrite(filename,heading10,Tname, 'A31');
          xlswrite(filename,heading5,Tname, 'A32');
          xlswrite(filename,heading11,Tname, 'A35');
          
          xlswrite(filename,heading12,Tname,'F2');
          xlswrite(filename,heading13,Tname, 'F3');
          xlswrite(filename,heading14,Tname, 'F4');
          xlswrite(filename,heading15,Tname, 'F5');
          xlswrite(filename,heading16,Tname, 'F6');
          xlswrite(filename,heading17,Tname, 'F7');
          xlswrite(filename,heading18,Tname, 'J1');
      
          xlswrite(filename,SLDJR_ARmom,Tname, 'B24:CX26');
          xlswrite(filename,SLDJR_ARmom_aver,Tname, 'B27:CX27');
          xlswrite(filename,SLDJR_KRmom,Tname, 'B28:CX30');
          xlswrite(filename,SLDJR_KRmom_aver,Tname,'B31:CX31');
          xlswrite(filename,SLDJR_HRmom,Tname, 'B32:CX34');
          xlswrite(filename,SLDJR_HRmom_aver,Tname, 'B35:CX35');
          
          
          
          xlswrite(filename,SLDJR_ARstiff,Tname, 'G5');
          xlswrite(filename,SLDJR_ARstiff_aver,Tname, 'J5');
          xlswrite(filename,SLDJR_KRstiff,Tname, 'G6');
          xlswrite(filename,SLDJR_KRstiff_aver,Tname, 'J6');
          xlswrite(filename,SLDJR_HRstiff,Tname, 'G7');
          xlswrite(filename,SLDJR_HRstiff_aver,Tname, 'J7');
         
         %COM Stiffness
          xlswrite(filename,heading19,Tname, 'L1');
          xlswrite(filename,COM_Stiff_SLDJR,Tname, 'L2');
          xlswrite(filename,heading18,Tname, 'O1');
          COM_Stiff_DJ_aver=mean(COM_Stiff_SLDJR);
          xlswrite(filename,COM_Stiff_DJ_aver,Tname, 'O2')
          
          
          
      end
      
  elseif contains(Tname, 'SLDJLeft')==1    
      AnkleLxSLDJL=interp1(1:numel(AnkleLxS),AnkleLxS, linspace(1, numel(AnkleLxS), 101))*1/1000;
      KneeLxSLDJL=interp1(1:numel(KneeLxS),KneeLxS, linspace(1, numel(KneeLxS), 101))*1/1000;
      HipLxSLDJL=interp1(1:numel(HipLxS),HipLxS, linspace(1, numel(HipLxS), 101))*1/1000;
      
      
      SLDJL_ALmom=[SLDJL_ALmom; AnkleLxSLDJL];
      
      SLDJL_KLmom=[SLDJL_KLmom; KneeLxSLDJL]; 
       
      SLDJL_HLmom=[SLDJL_HLmom; HipLxSLDJL];
      
      %*****************SLDJLeft Joint Stiffness***************       
        %Ankle Joint Stiffness
        %Left
       LAnkleAnglexIC=LAnkleAnglex(FS);
       LAnkleAnglexPeakFlex = LAnkleAnglexS(COMv0);
       Delta_LAnkleAngle=LAnkleAnglexPeakFlex-LAnkleAnglexIC
       Delta_AnkleLxS=AnkleLxS(COMv0)*1/1000;
       LAnkleJointStiff=Delta_AnkleLxS/Delta_LAnkleAngle;
        %Knee Joint Stiffness
        %Left
       LKneeAnglexIC=LKneeAnglex(FS);
       LKneeAnglexPeakFlex = LKneeAnglexS(COMv0);
       Delta_LKneeAngle=LKneeAnglexPeakFlex-LKneeAnglexIC;
       Delta_KneeLxS=KneeLxS(COMv0)*1/1000;
       LKneeJointStiff=Delta_KneeLxS/Delta_LKneeAngle;
        
        %Hip Joint Stiffness
        %Left
       LHipAnglexIC=LHipAnglex(FS);
       LHipAnglexPeakFlex =LHipAnglexS(COMv0);
       Delta_LHipAngle=LHipAnglexPeakFlex-LHipAnglexIC;
       Delta_HipLxS=HipLxS(COMv0)*1/1000;
       LHipJointStiff=Delta_HipLxS/Delta_LHipAngle;
       
       SLDJL_ALstiff=[SLDJL_ALstiff,LAnkleJointStiff];
      
       SLDJL_KLstiff=[SLDJL_KLstiff,LKneeJointStiff];
     
       SLDJL_HLstiff=[SLDJL_HLstiff,LHipJointStiff];
       
       %Calculating COM Stiffness
       COM_Stiff_SLDJL=[COM_Stiff_SLDJL,COM_Stiff];
        
      if contains(Tname,'SLDJLeft3')==1
          SLDJL_ALmom_aver=mean(SLDJL_ALmom);
          SLDJL_KLmom_aver=mean(SLDJL_KLmom);
          SLDJL_HLmom_aver=mean(SLDJL_HLmom);
          
          SLDJL_ALstiff_aver= mean(SLDJL_ALstiff);        
          SLDJL_KLstiff_aver= mean(SLDJL_KLstiff);
          SLDJL_HLstiff_aver= mean(SLDJL_HLstiff);
          
          
          xlswrite(filename,heading,Tname,'A12');
          xlswrite(filename,heading6,Tname,'A15');
          xlswrite(filename,heading1,Tname, 'A16');
          xlswrite(filename,heading7,Tname, 'A19');
          xlswrite(filename,heading2,Tname, 'A20');
          xlswrite(filename,heading8,Tname, 'A23');
          
          xlswrite(filename,heading12,Tname,'F2');
          xlswrite(filename,heading13,Tname, 'F3');
          xlswrite(filename,heading14,Tname, 'F4');
          xlswrite(filename,heading18,Tname, 'J1');
      
          xlswrite(filename,SLDJL_ALmom,Tname, 'B12:CX14');
          xlswrite(filename,SLDJL_ALmom_aver,Tname,'B15:CX15');
          xlswrite(filename,SLDJL_KLmom,Tname, 'B16:CX18');
          xlswrite(filename,SLDJL_KLmom_aver,Tname, 'B19:CX19');
          xlswrite(filename,SLDJL_HLmom,Tname, 'B20:CX22');
          xlswrite(filename,SLDJL_HLmom_aver,Tname, 'B23:CX23');
          
          xlswrite(filename,SLDJL_ALstiff,Tname, 'G2');
          xlswrite(filename,SLDJL_ALstiff_aver,Tname, 'J2');
          xlswrite(filename,SLDJL_KLstiff,Tname, 'G3');
          xlswrite(filename,SLDJL_KLstiff_aver,Tname, 'J3');
          xlswrite(filename,SLDJL_HLstiff,Tname, 'G4');
          xlswrite(filename,SLDJL_HLstiff_aver,Tname, 'J4');
         
         %COM Stiffness
          xlswrite(filename,heading19,Tname, 'L1');
          xlswrite(filename,COM_Stiff_SLDJL,Tname, 'L2');
          xlswrite(filename,heading18,Tname, 'O1');
          COM_Stiff_DJ_aver=mean(COM_Stiff_SLDJL);
          xlswrite(filename,COM_Stiff_DJ_aver,Tname, 'O2')
      end
      
xls='.xls';  
filename1=strcat(filename,xls);
  elseif contains(Tname, 'SLRJTRight')==1
      %Taking 5th rebound jump for moment analysis
      
      
      [num,txt]=xlsread(filename1 ,Tname,'B7:F7');
      txt=erase(txt,'Jump');
      txt=str2double(txt);
      for n=1:length(txt)
         AnkleRx_five=AnkleRx(FS(txt(n)):FO(txt(n)+1));
         AnkleRx_RJTR=interp1(1:numel(AnkleRx_five),AnkleRx_five, linspace(1, numel(AnkleRx_five), 101))*1/1000;
         AnkleRx_SLRJT=[AnkleRx_SLRJT;AnkleRx_RJTR];
         
         KneeRx_five=KneeRx(FS(txt(n)):FO(txt(n)+1));
         KneeRx_RJTR=interp1(1:numel(KneeRx_five),KneeRx_five, linspace(1, numel(KneeRx_five), 101))*1/1000;
         KneeRx_SLRJT=[KneeRx_SLRJT;KneeRx_RJTR];
         
         HipRx_five=HipRx(FS(txt(n)):FO(txt(n)+1));
         HipRx_RJTR=interp1(1:numel(HipRx_five),HipRx_five, linspace(1, numel(HipRx_five), 101))*1/1000;
         HipRx_SLRJT=[HipRx_SLRJT;HipRx_RJTR];
      end
      mean_AnkleRx_SLRJT=mean(AnkleRx_SLRJT);
      mean_KneeRx_SLRJT=mean(KneeRx_SLRJT);
      mean_HipRx_SLRJT=mean(HipRx_SLRJT);
      
      
      
      
      
      xlswrite(filename,heading3,Tname, 'A24');
      xlswrite(filename,heading4,Tname, 'A27');
      xlswrite(filename,heading5,Tname, 'A30');
      
      xlswrite(filename,mean_AnkleRx_SLRJT,Tname,'B24:CX24');
      xlswrite(filename,mean_KneeRx_SLRJT,Tname,'B27:CX27');
      xlswrite(filename,mean_HipRx_SLRJT,Tname,'B30:CX30');
      
      
      %*****************SLRJTRight Joint Stiffness***************       
        %Ankle Joint Stiffness
        %Right
      [num,txt]=xlsread(filename1 ,Tname,'B7:F7');
      txt=erase(txt,'Jump');
      txt=str2double(txt);
      for n=1:length(txt) 
       RAnkleAnglexIC=RAnkleAnglex(FS(txt(n)));
       AnkleRxS=AnkleRx(FS(txt(n)):FO(txt(n)+1));
       RAnkleAnglexS=RAnkleAnglex((FS(txt(n)):FO(txt(n)+1)));
       RAnkleAnglexPeakFlex = RAnkleAnglexS(COMv0_RJT(txt(n)));
       Delta_RAnkleAngle=RAnkleAnglexPeakFlex-RAnkleAnglexIC
       Delta_AnkleRxS=AnkleRxS(COMv0_RJT(txt(n)))*1/1000;
       RAnkleJointStiff=Delta_AnkleRxS/Delta_RAnkleAngle;
       RAnkleJointStiff_SLRJT=[RAnkleJointStiff_SLRJT,RAnkleJointStiff];
        
       
        %Knee Joint Stiffness
        %Right
       RKneeAnglexIC=RKneeAnglex(FS(txt(n)));
       KneeRxS=KneeRx(FS(txt(n)):FO(txt(n)+1));
       RKneeAnglexS=RKneeAnglex(FS(txt(n)):FO(txt(n)+1));
       RKneeAnglexPeakFlex =  RKneeAnglexS(COMv0_RJT(txt(n)));
       Delta_RKneeAngle=RKneeAnglexPeakFlex-RKneeAnglexIC
       Delta_KneeRxS=KneeRxS(COMv0_RJT(txt(n)))*1/1000;
       RKneeJointStiff=Delta_KneeRxS/Delta_RKneeAngle;
       RKneeJointStiff_SLRJT=[RKneeJointStiff_SLRJT,RKneeJointStiff];
        
       
        %Hip Joint Stiffness
        %Right
       RHipAnglexIC=RHipAnglex(FS(txt(n)));
       HipRxS=HipRx(FS(txt(n)):FO(txt(n)+1));
       RHipAnglexS=RHipAnglex((FS(txt(n)):FO(txt(n)+1)));
       RHipAnglexPeakFlex=RHipAnglexS(COMv0_RJT(txt(n)));
       Delta_RHipAngle=RHipAnglexPeakFlex-RHipAnglexIC;
       Delta_HipRxS=HipRxS(COMv0_RJT(txt(n)))*1/1000;
       RHipJointStiff=Delta_HipRxS/Delta_RHipAngle;
       RHipJointStiff_SLRJT=[RHipJointStiff_SLRJT,RHipJointStiff];
      end
      
      mean_RAnkleJointStiff_SLRJT= mean(RAnkleJointStiff_SLRJT);
      mean_RKneeJointStiff_SLRJT= mean(RKneeJointStiff_SLRJT);
      mean_RHipJointStiff_SLRJT= mean(RHipJointStiff_SLRJT);
      
      
       xlswrite(filename,heading15,Tname, 'N5'); 
       xlswrite(filename,heading16,Tname, 'N6');
       xlswrite(filename,heading17,Tname, 'N7');
      
       xlswrite(filename,mean_RAnkleJointStiff_SLRJT,Tname,'O5');
       xlswrite(filename,mean_RKneeJointStiff_SLRJT,Tname,'O6');
       xlswrite(filename,mean_RHipJointStiff_SLRJT,Tname,'O7');
       
       xlswrite(filename,heading19,Tname,'U1');
       xlswrite(filename,COM_Stiff_RJT(5),Tname,'U2');
       
  elseif contains(Tname, 'SLRJTLeft')==1
      %Taking 5th rebound jump for analysis
      [num,txt]=xlsread(filename1 ,Tname,'B7:F7');
      txt=erase(txt,'Jump');
      txt=str2double(txt);
      RKneeAnglexS=RKneeAnglex(FS(1):FO(2));
      LKneeAnglexS=LKneeAnglex(FS(1):FO(2));
      RAnkleAnglexS=RAnkleAnglex(FS(1):FO(2));
      LAnkleAnglexS=LAnkleAnglex(FS(1):FO(2));
      RHipAnglexS = RHipAnglex(FS(1):FO(2));
      LHipAnglexS=LHipAnglex(FS(1):FO(2));
 
      for n=1:length(txt)
         AnkleLx_five=AnkleLx(FS(txt(n)):FO(txt(n)+1));
         AnkleLx_RJTL=interp1(1:numel(AnkleLx_five),AnkleLx_five, linspace(1, numel(AnkleLx_five), 101))*1/1000;
         AnkleLx_SLRJT=[AnkleLx_SLRJT;AnkleLx_RJTL];
         
         KneeLx_five=KneeLx(FS(txt(n)):FO(txt(n)+1));
         KneeLx_RJTL=interp1(1:numel(KneeLx_five),KneeLx_five, linspace(1, numel(KneeLx_five), 101))*1/1000;
         KneeLx_SLRJT=[KneeLx_SLRJT;KneeLx_RJTL];
         
         HipLx_five=HipLx(FS(txt(n)):FO(txt(n)+1));
         HipLx_RJTL=interp1(1:numel(HipLx_five),HipLx_five, linspace(1, numel(HipLx_five), 101))*1/1000;
         HipLx_SLRJT=[HipLx_SLRJT;HipLx_RJTL];
      end
      mean_AnkleLx_SLRJT=mean(AnkleLx_SLRJT);
      mean_KneeLx_SLRJT=mean(KneeLx_SLRJT);
      mean_HipLx_SLRJT=mean(HipLx_SLRJT);
      
      xlswrite(filename,heading,Tname,'A12');
      xlswrite(filename,heading1,Tname, 'A15');
      xlswrite(filename,heading2,Tname, 'A18');
      
      xlswrite(filename,mean_AnkleLx_SLRJT,Tname,'B12:CX12');
      xlswrite(filename,mean_KneeLx_SLRJT,Tname,'B15:CX15');
      xlswrite(filename,mean_HipLx_SLRJT,Tname,'B18:CX18');
      
      %*****************SLRJTLeft Joint Stiffness***************       
         %Ankle Joint Stiffness
      [num,txt]=xlsread(filename1 ,Tname,'B7:F7');
      txt=erase(txt,'Jump');
      txt=str2double(txt);
      for n=1:length(txt) 
       LAnkleAnglexIC=LAnkleAnglex(FS(txt(n)));
       AnkleLxS=AnkleLx(FS(txt(n)):FO(txt(n)+1));
       LAnkleAnglexS=LAnkleAnglex(FS(txt(n)):FO(txt(n)+1));
       LAnkleAnglexPeakFlex = LAnkleAnglexS(COMv0_RJT(txt(n)));
       Delta_LAnkleAngle=LAnkleAnglexPeakFlex-LAnkleAnglexIC
       Delta_AnkleLxS=AnkleLxS(COMv0_RJT(txt(n)))*1/1000;
       LAnkleJointStiff=Delta_AnkleLxS/Delta_LAnkleAngle;
       LAnkleJointStiff_SLRJT=[LAnkleJointStiff_SLRJT,LAnkleJointStiff];
       
       %Knee Joint Stiffness
        
       LKneeAnglexIC=LKneeAnglex(FS(txt(n)));
       KneeLxS=KneeLx(FS(txt(n)):FO(txt(n)+1));
       LKneeAnglexS=LKneeAnglex((FS(txt(n)):FO(txt(n)+1)));
       LKneeAnglexPeakFlex = LKneeAnglexS(COMv0_RJT(txt(n)));
       Delta_LKneeAngle=LKneeAnglexPeakFlex-LKneeAnglexIC;
       Delta_KneeLxS=KneeLxS(COMv0_RJT(txt(n)))*1/1000;
       LKneeJointStiff=Delta_KneeLxS/Delta_LKneeAngle;
       LKneeJointStiff_SLRJT= [LKneeJointStiff_SLRJT,LKneeJointStiff];
       
       %Hip Joint Stiffness
        
       LHipAnglexIC=LHipAnglex(FS(txt(n)));
       HipLxS=HipLx(FS(txt(n)):FO(txt(n)+1));
       LHipAnglexS=LHipAnglex((FS(txt(n)):FO(txt(n)+1)));
       LHipAnglexPeakFlex =LHipAnglexS(COMv0_RJT(txt(n)));
       Delta_LHipAngle=LHipAnglexPeakFlex-LHipAnglexIC;
       Delta_HipLxS=HipLxS(COMv0_RJT(txt(n)))*1/1000;
       LHipJointStiff=Delta_HipLxS/Delta_LHipAngle;
       LHipJointStiff_SLRJT=[LHipJointStiff_SLRJT,LHipJointStiff];
      end 
       
      mean_LAnkleJointStiff_SLRJT=mean(LAnkleJointStiff_SLRJT);
      mean_LKneeJointStiff_SLRJT=mean(LKneeJointStiff_SLRJT);
      mean_LHipJointStiff_SLRJT=mean(LHipJointStiff_SLRJT);
      
       xlswrite(filename,heading12,Tname, 'N2'); 
       xlswrite(filename,heading13,Tname, 'N3');
       xlswrite(filename,heading14,Tname, 'N4');
      
       xlswrite(filename,mean_LAnkleJointStiff_SLRJT,Tname,'O2');
       xlswrite(filename,mean_LKneeJointStiff_SLRJT,Tname,'O3');
       xlswrite(filename,mean_LHipJointStiff_SLRJT,Tname,'O4');
       
       xlswrite(filename,heading19,Tname,'U1');
       xlswrite(filename,COM_Stiff_RJT(5),Tname,'U2');
      
  else
      %Taking 5th rebound jump from DLRJT
      xls='.xls';  
      filename1=strcat(filename,xls);
      [num,txt]=xlsread(filename1 ,'RJT','B7:F7');
      txt=erase(txt,'Jump');
      txt=str2double(txt);
      for n=1:length(txt)
        AnkleRx_five=AnkleRx(FS(txt(n)):FO(txt(n)+1));
        AnkleRx_RJTR=interp1(1:numel(AnkleRx_five),AnkleRx_five, linspace(1, numel(AnkleRx_five), 101))*1/1000;
        AnkleRx_DLRJT=[AnkleRx_DLRJT;AnkleRx_RJTR];
         
        KneeRx_five=KneeRx(FS(txt(n)):FO(txt(n)+1));
        KneeRx_RJTR=interp1(1:numel(KneeRx_five),KneeRx_five, linspace(1, numel(KneeRx_five), 101))*1/1000;
        KneeRx_DLRJT=[KneeRx_DLRJT;KneeRx_RJTR];
         
        HipRx_five=HipRx(FS(txt(n)):FO(txt(n)+1));
        HipRx_RJTR=interp1(1:numel(HipRx_five),HipRx_five, linspace(1, numel(HipRx_five), 101))*1/1000;
        HipRx_DLRJT=[HipRx_DLRJT;HipRx_RJTR];
      
        AnkleLx_five=AnkleLx(FS(txt(n)):FO(txt(n)+1));
        AnkleLx_RJTL=interp1(1:numel(AnkleLx_five),AnkleLx_five, linspace(1, numel(AnkleLx_five), 101))*1/1000;
        AnkleLx_DLRJT=[AnkleLx_DLRJT;AnkleLx_RJTL];
         
        KneeLx_five=KneeLx(FS(txt(n)):FO(txt(n)+1));
        KneeLx_RJTL=interp1(1:numel(KneeLx_five),KneeLx_five, linspace(1, numel(KneeLx_five), 101))*1/1000;
        KneeLx_DLRJT=[KneeLx_DLRJT;KneeLx_RJTL];
         
        HipLx_five=HipLx(FS(txt(n)):FO(txt(n)+1));
        HipLx_RJTL=interp1(1:numel(HipLx_five),HipLx_five, linspace(1, numel(HipLx_five), 101))*1/1000;
        HipLx_DLRJT=[HipLx_DLRJT;HipLx_RJTL];
      end
      
      mean_AnkleRx_DLRJT=mean(AnkleRx_DLRJT);
      mean_KneeRx_DLRJT=mean(KneeRx_DLRJT);
      mean_HipRx_DLRJT=mean(HipRx_DLRJT);
      
      mean_AnkleLx_DLRJT=mean(AnkleLx_DLRJT);
      mean_KneeLx_DLRJT=mean(KneeLx_DLRJT);
      mean_HipLx_DLRJT=mean(HipLx_DLRJT);
      
      xlswrite(filename,heading,Tname,'A12');
      xlswrite(filename,heading1,Tname, 'A15');
      xlswrite(filename,heading2,Tname, 'A18');
      xlswrite(filename,heading3,Tname, 'A24');
      xlswrite(filename,heading4,Tname, 'A27');
      xlswrite(filename,heading5,Tname, 'A30');
      
      xlswrite(filename,mean_AnkleLx_DLRJT,Tname,'B12:CX12');
      xlswrite(filename,mean_KneeLx_DLRJT,Tname,'B15:CX15');
      xlswrite(filename,mean_HipLx_DLRJT,Tname,'B18:CX18');
      xlswrite(filename,mean_AnkleRx_DLRJT,Tname,'B24:CX24');
      xlswrite(filename,mean_KneeRx_DLRJT,Tname,'B27:CX27');
      xlswrite(filename,mean_HipRx_DLRJT,Tname,'B30:CX30');
      
      %*****************RJT Joint Stiffness***************       
         %Ankle Joint Stiffness
         %Left
      [num,txt]=xlsread(filename1 ,Tname,'B7:F7');
      txt=erase(txt,'Jump');
      txt=str2double(txt);
      for n=1:length(txt) 
       LAnkleAnglexIC=LAnkleAnglex(FS(txt(n)));
       AnkleLxS=AnkleLx(FS(txt(n)):FO(txt(n)+1));
       LAnkleAnglexS=LAnkleAnglex((FS(txt(n)):FO(txt(n)+1)));
       LAnkleAnglexPeakFlex = LAnkleAnglexS(COMv0_RJT(txt(n)));
       Delta_LAnkleAngle=LAnkleAnglexPeakFlex-LAnkleAnglexIC
       Delta_AnkleLxS=AnkleLxS(COMv0_RJT(txt(n)))*1/1000;
       LAnkleJointStiff=Delta_AnkleLxS/Delta_LAnkleAngle;
       LAnkleJointStiffness_DLRJT=[LAnkleJointStiffness_DLRJT,LAnkleJointStiff];
        
       
        %Knee Joint Stiffness
        %left
       LKneeAnglexIC=LKneeAnglex(FS(txt(n)));
       KneeLxS=KneeLx(FS(txt(n)):FO(txt(n)+1));
       LKneeAnglexS=LKneeAnglex(FS(txt(n)):FO(txt(n)+1));
       LKneeAnglexPeakFlex = LKneeAnglexS(COMv0_RJT(txt(n)));
       Delta_LKneeAngle=LKneeAnglexPeakFlex-LKneeAnglexIC
       Delta_KneeLxS=KneeLxS(COMv0_RJT(txt(n)))*1/1000;
       LKneeJointStiff=Delta_KneeLxS/Delta_LKneeAngle;
       LKneeJointStiffness_DLRJT=[LKneeJointStiffness_DLRJT,LKneeJointStiff];
        
       
        %Hip Joint Stiffness
        %Left
       LHipAnglexIC=LHipAnglex(FS(txt(n)));
       HipLxS=HipLx(FS(txt(n)):FO(txt(n)+1));
       LHipAnglexS=LHipAnglex(FS(txt(n)):FO(txt(n)+1));
       LHipAnglexPeakFlex = LHipAnglexS(COMv0_RJT(txt(n)));
       Delta_LHipAngle=LHipAnglexPeakFlex-LHipAnglexIC
       Delta_HipLxS=HipLxS(COMv0_RJT(txt(n)))*1/1000;
       LHipJointStiff=Delta_HipLxS/Delta_LHipAngle;
       LHipJointStiffness_DLRJT=[LHipJointStiffness_DLRJT,LHipJointStiff];
       
        %Ankle Joint Stiffness
        %Right
       RAnkleAnglexIC=RAnkleAnglex(FS(txt(n)));
       AnkleRxS=AnkleRx(FS(txt(n)):FO(txt(n)+1));
       RAnkleAnglexS=RAnkleAnglex(FS(txt(n)):FO(txt(n)+1));
       RAnkleAnglexPeakFlex =RAnkleAnglexS(COMv0_RJT(txt(n)));
       Delta_RAnkleAngle=RAnkleAnglexPeakFlex-RAnkleAnglexIC
       Delta_AnkleRxS=AnkleRxS(COMv0_RJT(txt(n)))*1/1000;
       RAnkleJointStiff=Delta_AnkleRxS/Delta_RAnkleAngle;
       RAnkleJointStiffness_DLRJT=[RAnkleJointStiffness_DLRJT,RAnkleJointStiff];
        
       
        %Knee Joint Stiffness
        %Right
       RKneeAnglexIC=RKneeAnglex(FS(txt(n)));
       KneeRxS=KneeRx(FS(txt(n)):FO(txt(n)+1));
       RKneeAnglexS=RKneeAnglex(FS(txt(n)):FO(txt(n)+1));
       RKneeAnglexPeakFlex=RKneeAnglexS(COMv0_RJT(txt(n)));
       Delta_RKneeAngle=RKneeAnglexPeakFlex-RKneeAnglexIC
       Delta_KneeRxS=KneeRxS(COMv0_RJT(txt(n)))*1/1000;
       RKneeJointStiff=Delta_KneeRxS/Delta_RKneeAngle;
       RKneeJointStiffness_DLRJT=[RKneeJointStiffness_DLRJT,RKneeJointStiff];
        
       
        %Hip Joint Stiffness
        %Right
       RHipAnglexIC=RHipAnglex(FS(txt(n)));
       HipRxS=HipRx(FS(txt(n)):FO(txt(n)+1));
       RHipAnglexS=RHipAnglex(FS(txt(n)):FO(txt(n)+1));
       RHipAnglexPeakFlex =RHipAnglexS(COMv0_RJT(txt(n)));
       Delta_RHipAngle=RHipAnglexPeakFlex-RHipAnglexIC
       Delta_HipRxS=HipRxS(COMv0_RJT(txt(n)))*1/1000;
       RHipJointStiff=Delta_HipRxS/Delta_RHipAngle;
       RHipJointStiffness_DLRJT=[RHipJointStiffness_DLRJT,RHipJointStiff];
      end   
      
      mean_RAnkleJointStiffness_DLRJT=mean(RAnkleJointStiffness_DLRJT);
      mean_RKneeJointStiffness_DLRJT=mean(RKneeJointStiffness_DLRJT);
      mean_RHipJointStiffness_DLRJT=mean(RHipJointStiffness_DLRJT);
      mean_LAnkleJointStiffness_DLRJT=mean(LAnkleJointStiffness_DLRJT);
      mean_LKneeJointStiffness_DLRJT=mean(LKneeJointStiffness_DLRJT);
      mean_LHipJointStiffness_DLRJT=mean(LHipJointStiffness_DLRJT);
     
       
       xlswrite(filename,heading12,Tname, 'N2'); 
       xlswrite(filename,heading13,Tname, 'N3');
       xlswrite(filename,heading14,Tname, 'N4');
       xlswrite(filename,heading15,Tname, 'N5'); 
       xlswrite(filename,heading16,Tname, 'N6');
       xlswrite(filename,heading17,Tname, 'N7');
      
      
       xlswrite(filename,mean_LAnkleJointStiffness_DLRJT,Tname,'O2');
       xlswrite(filename,mean_LKneeJointStiffness_DLRJT,Tname,'O3');
       xlswrite(filename,mean_LHipJointStiffness_DLRJT,Tname,'O4');
       xlswrite(filename,mean_RAnkleJointStiffness_DLRJT,Tname,'O5');
       xlswrite(filename,mean_RKneeJointStiffness_DLRJT,Tname,'O6');
       xlswrite(filename,mean_RHipJointStiffness_DLRJT,Tname,'O7');
       
       xlswrite(filename,heading19,Tname,'U1');
       xlswrite(filename,COM_Stiff_RJT(5),Tname,'U2');
  end
  
  
end
    
if injside==1
   LSI=(RJT_Data_L./RJT_Data_R)*100;
   LSI1= round(LSI,2);
   xlswrite(filename,LSI1,'SLRJTLeft','I8')
   LSI2= (AVER_L./AVER_R)*100;
   LSI3= round(LSI2,2);
   xlswrite(filename,LSI3,'SLDJLeft3','E5')
   LSI_header={'LSI%'};
   xlswrite(filename,LSI_header,'SLDJLeft3','E4')
else
   LSI=(RJT_Data_R./RJT_Data_L)*100;
   LSI1= round(LSI,2);
   xlswrite(filename,LSI1,'SLRJTRight','I8')
   LSI2= (AVER_R./AVER_L)*100;
   LSI3= round(LSI2,2);
   xlswrite(filename,LSI3,'SLDJRight3','E5')
   LSI_header={'LSI%'};
   xlswrite(filename,LSI_header,'SLDJRight3','E4')
end


    
    
   
    
 
    
   
        
   
   
   
    
    %stop loop if only one file is selected
    if ischar(Fname)
        return
    end   
    
    
     
    
 


