function [dataStreamOut, bitTime] = manchesterDecodeFloat2(dataStreamIn, rawTime, resyncThreshold)
%CONSTELLATION - convert to bits from raw manchester bits
%resyncThreshold = 1;

%Idea - this threshold could be made dynamic and could be a percentage of
%a moving average over the last X points. This would let the algorithm
%follow gain->amplitude variations caused by interference from nearby signals.

idx2=1;
clockmod = 0;
idxerr=1;
%clear errx erry dataStreamOut BitTime;
dataStreamOut(1) = uint8(0);
prevSample = 0;
currentSample = 0;
evenOddCounter = 0;

for idx = 1:numel(dataStreamIn)
    prevPrevSample = prevSample;
    prevSample = currentSample;
    currentSample = dataStreamIn(idx);
    %If not a bit boundary, see if it should be and we're out of sync
    %But only resync on strong bits
    if(mod(evenOddCounter,2) ~= clockmod)
        if(sign(prevPrevSample) == sign(prevSample))
            errx(idxerr)=idx2;
            erry(idxerr)=prevSample;
            idxerr=idxerr+1;
            if(abs(prevPrevSample) > resyncThreshold && abs(prevSample) > resyncThreshold)
                clockmod = mod(evenOddCounter,2); %only resync if we have confidence in BOTH bit decisions
            end
            
        end
    end
    
    %check for bit boundary, and make decision using the strongest of the
    %two bits.
    if(mod(evenOddCounter,2) == clockmod)
        if(abs(prevSample) > abs(currentSample)) %use the strongest symbol to determine bit
            if(prevSample > 0)
                currentBit = '1';
            else
                currentBit = '0';
            end
        else
            if(currentSample > 0)
                currentBit = '0';
            else
                currentBit = '1';
            end
        end
        dataStreamOut(idx2) = currentBit;
        bitTime(idx2)=rawTime(idx);
        idx2 = idx2+1;
    end
    
    evenOddCounter = evenOddCounter+1;
end
fprintf([num2str(idxerr) ' errors\n']);
end