function [FS, FO] = getevents(forcedata, BW)
%this function is designed to use the vertical
%force date to calculate the footstrike and footoff
%events of up to 10 repeated jumps
%Author - Yunus Kaghembe

%find how many hops participants have done
[peaks, loc] = findpeaks(forcedata, 'MinPeakHeight', BW + 200, 'MinPeakDistance', 100);

%Classifying start of jumps/hops
Start = find(forcedata < BW - 40, 1, 'first');
FootOff = [];
FootStrike =[];
FS = [];
FO = [];
if length(loc) > 11
    for  z = 1:11
        if isempty(FootOff)
            FootOff = find(forcedata(Start:end,:) < 0, 1, 'first') + Start; %First FootOff
            FO(:,z)=FootOff;
        elseif ~isempty(FootOff)
            FootOff = find(forcedata(FootStrike:end,:) < 0, 1, 'first') + FootStrike;
            FO(:,z)=FootOff;
        end
        FootStrike = find(forcedata(FootOff:end, :) >= 20, 1, 'first') + FootOff - 1;
        FS(:,z) = FootStrike;
    end
elseif length(loc) <=2
    for z=1:length(loc)-1
        if isempty(FootStrike)
            FootStrike= find(forcedata(Start:loc, :)>= 20,1, 'first') + Start
            FS(:,z)=FootStrike;
        end
        FootOff= find(forcedata(FootStrike:end, :) < 0, 1, 'first') + FootStrike
        FO(:,z)=FootOff;
    end   
elseif numel(loc)<=3
    for z=1:length(loc)-2
        if isempty(FootStrike)
            FootStrike= find(forcedata(Start:loc, :)>= 20,1, 'first') + Start
            FS(:,z)=FootStrike;
        end
        FootOff= find(forcedata(FootStrike:end, :) < 0, 1, 'first') + FootStrike
        FO(:,z)=FootOff;
    end   
    
elseif length(loc) <= 11 
    for  z = 1:length(loc)-1
        if isempty(FootOff)
            FootOff = find(forcedata(Start:end,:) < 0, 1, 'first') + Start; %First FootOff
            FO(:,z)=FootOff;
        elseif ~isempty(FootOff)
            FootOff = find(forcedata(FootStrike:end,:) < 0, 1, 'first') + FootStrike;
            FO(:,z)=FootOff;
        end
        FootStrike = find(forcedata(FootOff:end, :) >= 20, 1, 'first') + FootOff - 1;
        FS(:,z) = FootStrike;
    end

            
end
end

