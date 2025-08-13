function [ContactTime,JumpHeightDJ,JumpHeightv1,RSI] = RJTanalysis(FS,FO,CoMz)
    
%this fuctions used identified Footstrike (FS)  and  Footoff (FO)
%events to calculate contact time, jump height and reactive strength index..
%(RSI)
%Author - Yunus Kaghembe & Leigh J Ryan 

%Calculating Variables

a=double(FS);
b=double(FO);
if length(b) < 2
    ab = length(b)
elseif length(b) < 11
    ab = length(b)-1;
elseif length (b) >= 11
    ab = 10;
end
FT1=[];
for h = 1:ab
    %Frametime1= a-b;
    if length(b) < 2
        Frametime= b-a
        FT1=Frametime
    else    
    Frametime=b(:,h+1)-a(:,h);
    FT1(:,h)=Frametime;
    end
end
ContactTime = round((FT1/200),3); %Contact Time

%Calculating Flight time
%if length(a) < 11
    %xd = length(a);
    %c=a(2:xd);
    %d=b(2:xd);
%elseif length(a) >= 11
    %c=a(2:11);
    %d=b(2:11);
%end
%Frametime2=c-d;
%FlightTime=round((Frametime2/200),2);
%vv= (9.81*FlightTime)./2;
%JumpHeight= round(((vv.^2)/(2*9.81)),2);

%end%
%BW=mean(VerticalForce(1:100));

%BM=BW/9.81;
%I=zeros(size(FBW));
%dt=1/200;
%for n=1:length(VerticalForce)-1
    %I(n+1)=((FBW(n)+FBW(n+1))/2)*dt+I(n);
%end 
%A=[]
%JumpHeight_all=[]
%vel=I/BM;
%for i=1:(length(FO)-1)
    %Takeoff_vel=vel(FO(i+1));
    %A=[A,Takeoff_vel];
   
    %JumpHeight=round(((A.^2)/(2*9.81)),3);
    

    






   

%Calculating jump heigth for RJT AND DJ


COMv=(diff(CoMz)*Fs)/1000;
JumpHeightv1=[];
JumpHeightDJ=[];
if length(FO)<2
    JumpHeightDJ=(((COMv(FO)).^2)/(2*9.81));
    t=1:1:length(COMz)-1
    JumpHeightDJ=[JumpHeightDJ];
else    
    for n=1:length(FO)-1
    JumpHeightv=((COMv(FO(n)).^2)/(2*9.81));
    JumpHeightv1=[JumpHeightv1,JumpHeightv];
    t=1:1:length(COMz)-1;
 
    end
end
if length(FO)<2
    RSI = round((JumpHeightDJ./ContactTime),3);
else    
    RSI = round((JumpHeightv1./ContactTime),3);
end
plot(t,COMv)

if length(FO)<2
    plot(t(FO),COMv(FO),'r*')
else    
    plot(t(FO(1:length(FO)-1)),COMv(FO(1:length(FO)-1)),'r*')
end    
 
%pause;
%close
end


