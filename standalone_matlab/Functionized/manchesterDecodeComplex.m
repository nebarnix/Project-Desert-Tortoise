function [dataStreamOut] = manchesterDecode(dataStreamIn,resyncThreshold)
%CONSTELLATION - convert to bits from raw manchester bits
%resyncThreshold = 1;

%Idea - this threshold could be made dynamic and could be a percentage of
%a moving average over the last X points. This would let the algorithm
%follow gain->amplitude variations caused by interference from nearby signals. 

idx2=1;
clockmod = 0;
idxerr=1;
clear errx erry dataStreamOut BitTime;
dataStreamOut(1) = uint8(0);

for idx = 2:numel(dataStreamIn)-1    
    %If not a bit boundary, see if it should be and we're out of sync
    %But only resync on strong bits
    if(mod(idx,2) ~= clockmod)    
        if(sign(imag(dataStreamIn(idx-1))) == sign(imag(dataStreamIn(idx))))
            errx(idxerr)=idx2;
            erry(idxerr)=imag(dataStreamIn(idx));
            idxerr=idxerr+1;   
            if(abs(imag(dataStreamIn(idx-1))) > resyncThreshold && abs(imag(dataStreamIn(idx))) > resyncThreshold)                
                clockmod = mod(idx,2); %only resync if we have confidence in BOTH bit decisions
            end
            
        end        
    end
    
    %check for bit boundary, and make decision using the strongest of the
    %two bits. 
    if(mod(idx,2) == clockmod)
        if(abs(imag(dataStreamIn(idx))) > abs(imag(dataStreamIn(idx+1)))) %use the strongest symbol to determine bit
            if(imag(dataStreamIn(idx)) > 0)                
                dataStreamOut(idx2) = '0';
            else
                dataStreamOut(idx2) = '1';
            end
        else
            if(imag(dataStreamIn(idx+1)) > 0)                
                dataStreamOut(idx2) = '1';
            else
                dataStreamOut(idx2) = '0';
            end        
        end
        
        BitTime(idx2)=RawTime(idx);    
        idx2 = idx2+1;           
    end        
end