%M&M Clock Recovery Loop (interpolating version!)
function [dataStreamOut, dataStreamOutTime] = MMClockRecovery(dataStreamIn, Fs, baud, stepSpread)

%baud=8320*2-1;
%stepSpread=10;
stepSize = Fs/(baud);
stepMax = Fs/(baud-stepSpread);
stepMin = Fs/(baud+stepSpread);

%dataStreamOut = zeros(1,round(length(y)/stepSize));
%Ind  = zeros(1,round(length(y)/stepSize));

kp = 0.25;
%kp = .025;
Count = 1;
nextSample = 1;
sampleLast  = 1;
while nextSample < length(dataStreamIn)
    
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