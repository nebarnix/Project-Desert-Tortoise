%M&M Clock Recovery Loop (interpolating version!)
function [dataStreamOut, dataStreamOutTime] = UpsamplingGardenerClockRecovery(dataStreamIn, dataStreamInTime, FsIn, FsOut, baud, stepSpread, kp)

%FsOut=8320*15;
%FsOut=50e3;
fprintf('Interpolating Up...');
interpLevel = FsOut/FsIn;
%dataStreamIn = interp1(real(dataStreamIn),1:1/interpLevel:numel(real(dataStreamIn)),'linear');
%dataStreamIn = interp1(real(dataStreamIn),1:1/interpLevel:numel(real(dataStreamIn)),'spline');
dataStreamIn = interp1(real(dataStreamIn),1:1/interpLevel:numel(real(dataStreamIn)),'cubic');

%mmInDataStream = mmInDataStream(22976652:24565761);

Fs=FsOut;
Ts = 1/(Fs);
dataStreamInTime=0:Ts:Ts*(numel(dataStreamIn)-1);

fprintf('done\n');

%baud=8320*2-1;
%stepSpread=10;
stepSize = Fs/(baud);
%stepMax = Fs/(baud-stepSpread);
%stepMin = Fs/(baud+stepSpread);

%dataStreamOut = zeros(1,round(length(y)/stepSize));
%Ind  = zeros(1,round(length(y)/stepSize));

%kp = 0.25;
%kp = .025;
Count = 1;
nextSample = 1;
prevBit  = 1;
halfSample  = 1;
progress = 0;
idx = 1;
percentcomplete=0;
onePercent = (numel(dataStreamIn)/stepSize) / 100;
msgCount = 1;
fprintf('MM Clock Recovery:');
MSG = [num2str(percentcomplete) '%%'];
        msgCount = numel(MSG);
        fprintf(MSG);

while nextSample < length(dataStreamIn)
    
    progress = progress+1;
    if(progress >= 1*onePercent)        
        for progress=1:msgCount-1
           fprintf(char(8));
        end
        percentcomplete = percentcomplete + 1;
        MSG = [num2str(percentcomplete) '%%'];
        msgCount = numel(MSG);
        fprintf(MSG);        
        progress = 0;
    end
    
    %Stores Bit 
    currentBit  = dataStreamIn(round(nextSample));
    halfSample  = dataStreamIn(round(halfSample));
    dataStreamOutTime(Count) = dataStreamInTime(round(nextSample));
    dataStreamOut(Count) = currentBit;
    %Ind(Count)  = nextSample;
    Count   = Count + 1;
    
    %Calculates Error
    %Error = sign(real(prevBit))*real(currentBit) - sign(real(currentBit))*real(prevBit);
    Error = kp*(currentBit - prevBit).*(halfSample);
    %if error is LESS THAN 0 a timing advance is required (stepsize)
    %if error is GREATER THAN 0 a timing delay is required (stepsize)
   
    %Updates Step Size
    %stepSize = stepSize + kp*Error;
      
    %Phase = kp*Error;
    
    
    %Limits Step size
    if( Error > stepSpread)
        Error = stepSpread;
    end
    if( Error < -stepSpread )
        Error = -stepSpread;
    end
    
    errori(idx) = Error;
    %stepSizei(idx) = Error;
    nextSample = nextSample - Error;
    
    
    %Updates nextSample
    halfSample = nextSample + stepSize/2.0;
    nextSample = nextSample + stepSize;    
    prevBit = currentBit;
    idx= idx+1;
end
for progress=1:msgCount-1      
    fprintf(char(8));        
end
fprintf('100%%\n');
%subplot(2,1,1);
plot(errori);
%subplot(2,1,2);
%plot(stepSizei);
end