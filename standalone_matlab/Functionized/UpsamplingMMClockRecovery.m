%M&M Clock Recovery Loop (interpolating version!)
function [dataStreamOut, dataStreamOutTime] = UpsamplingMMClockRecovery(dataStreamIn, dataStreamInTime, FsIn, FsOut, baud, stepSpread, kp)

%FsOut=8320*15;
%FsOut=50e3;
fprintf('Interpolating Up...');
interpLevel = FsOut/FsIn;
dataStreamIn = interp1(real(dataStreamIn),1:1/interpLevel:numel(real(dataStreamIn)),'spline');
%mmInDataStream = mmInDataStream(22976652:24565761);

Fs=FsOut;
Ts = 1/(Fs);
dataStreamInTime=0:Ts:Ts*(numel(dataStreamIn)-1);

fprintf('done\n');

%baud=8320*2-1;
%stepSpread=10;
stepSize = Fs/(baud);
stepMax = Fs/(baud-stepSpread);
stepMin = Fs/(baud+stepSpread);

%dataStreamOut = zeros(1,round(length(y)/stepSize));
%Ind  = zeros(1,round(length(y)/stepSize));

%kp = 0.25;
%kp = .025;
Count = 1;
nextSample = 1;
sampleLast  = 1;

progress = 0;
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
    dataStreamOutTime(Count) = dataStreamInTime(round(nextSample));
    dataStreamOut(Count) = currentBit;
    %Ind(Count)  = nextSample;
    Count   = Count + 1;
    
    %Calculates Error
    Error = sign(real(sampleLast))*real(currentBit) - sign(real(currentBit))*real(sampleLast);
    
    %Updates Step Size
    stepSize = stepSize + kp*Error;
    
    %Limits Step size
    if( stepSize > stepMax )
        stepSize = stepMax;
    end
    if( stepSize < stepMin )
        stepSize = stepMin;
    end
    
    %Updates nextSample
    nextSample = nextSample + stepSize;
    sampleLast = currentBit;
    
end
for progress=1:msgCount-1      
    fprintf(char(8));        
end
fprintf('100%%\n');
end