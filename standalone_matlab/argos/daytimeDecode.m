%BITSTREAM - Print Major Frame UTC times and spacecraft ID
function [dayNum, spaceCraft, minorFrameID] = daytimeDecode(minorFrames, frameTime)
clear minorFrameID dayNum spacecraftID;
minorFrameID = bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6)); %pull out minor frame counter, which is annoyingly 9 bits long

%If you want plots, uncomment below
%figure(4);
%plot(frameTime(:,5),minorFrameID,'.'); %plot minor frame counter

spaceCraft(1) = 0;
dayNum(1) = 0;

idx = 1;

for frame=1:numel(minorFrameID)        
    spaceCraft(frame) = minorFrames(frame,3);    
    
    if(minorFrameID(frame) == 0)
        dayNum(idx)=bitshift(minorFrames(frame,8+1),1)+bitshift(bitor(minorFrames(frame,9+1),128),-7); %should be 241 for loproto2 recording        
        
        %check if the counter is more than the number of ms per day (bad data)
        if((bitshift(bitand(minorFrames(frame,9+1),bin2dec('111')),24) + bitshift(minorFrames(frame,10+1),16) + bitshift(minorFrames(frame,11+1),8) + minorFrames(frame,12+1)) < 86400000)
            dayMSeconds(idx) = (bitshift(bitand(minorFrames(frame,9+1),bin2dec('111')),24) + bitshift(minorFrames(frame,10+1),16) + bitshift(minorFrames(frame,11+1),8) + minorFrames(frame,12+1));
            fprintf([num2str(frameTime(frame,1)) ' Local Seconds is ' num2str(dayMSeconds(idx)/1000.0) ' Spacecraft Day Seconds' ]);
            
            hour = floor(dayMSeconds(idx)/(1000*60*60));
            minute = floor((dayMSeconds(idx)/(1000*60*60) - hour)*60);
            seconds = (((dayMSeconds(idx)/(1000*60*60) - hour)*60) - minute) *60;
            fprintf([' which is ' num2str(hour) ':' num2str(minute) ':' num2str(seconds)]);
            
            if idx > 1 && dayMSeconds(idx) <= dayMSeconds(idx-1)
                fprintf(' ...but this might be an error');
            else
                hour = floor((dayMSeconds(idx)-frameTime(frame,1)*1000)/(1000*60*60));
                minute = floor(((dayMSeconds(idx)-frameTime(frame,1)*1000)/(1000*60*60) - hour)*60);
                seconds = ((((dayMSeconds(idx)-frameTime(frame,1)*1000)/(1000*60*60) - hour)*60) - minute) *60;
                fprintf(['    t(0) => ' num2str(hour) ':' num2str(minute) ':' num2str(seconds)]);
            end
            fprintf('\n');
        else
            dayMSeconds(idx) = -1;
        end        
        idx = idx + 1;
        %fprintf([num2str(frameTime(frame,1)) ':' num2str(dec2bin(minorFrames(frame,9+1),8)) ' ' num2str(dec2bin(minorFrames(frame,10+1),8)) ' ' num2str(dec2bin(minorFrames(frame,11+1),8)) ' ' num2str(dec2bin(minorFrames(frame,12+1),8)) '\n']);
        %fprintf([num2str(day) ' ' num2str(dec2hex(minorFrames(frame,8+1))) ' ' num2str(dec2hex(minorFrames(frame,9+1))) '\n\n']);        
    end
end

if mode(spaceCraft) == 8
   fprintf(['Spacecraft: ' num2str(mode(spaceCraft)) '=>NOAA-15\n']);    
   fprintf(['Julean Day: ' num2str(mode(dayNum)) ' \n']);
elseif mode(spaceCraft) == 13
    fprintf(['Spacecraft: ' num2str(mode(spaceCraft)) '=>NOAA-18\n']);    
    fprintf(['Julean Day: ' num2str(mode(dayNum)) ' \n']);
elseif mode(spaceCraft) == 15
    fprintf(['Spacecraft: ' num2str(mode(spaceCraft)) '=>NOAA-19\n']);    
    fprintf(['Julean Day: ' num2str(mode(dayNum)) ' \n']);
else
    fprintf(['Spacecraft: ' num2str(mode(spaceCraft)) '=>A UFO!\n']);    
end
spaceCraft = mode(spaceCraft);
dayNum = mode(dayNum);
end