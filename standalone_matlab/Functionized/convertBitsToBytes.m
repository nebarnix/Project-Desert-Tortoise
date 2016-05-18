%convert ascii binary stream to actual binary at matched syncword locations
%Includes correct and inverted syncwords for periods of constellation reversal. 

function [minorFrames, frameTime] = convertBitsToBytes(dataStreamIn, bitTime, SyncWordIndex, SyncWordInvIndex)
clear minorFrames frameTime;

SyncWordAllIndex = sort(cat(2,SyncWordIndex,SyncWordInvIndex));

for frameIdx=1:numel(SyncWordAllIndex)-1
    %See if the frame is normal or inverted bits
    if isempty(find(SyncWordInvIndex == SyncWordAllIndex(frameIdx),1))
        for frameByteIdx=0:103 %minor frames are 103 bytes long
            byte=0;
            %Start of byte time
            frameTime(frameIdx,frameByteIdx+1)=bitTime(SyncWordAllIndex(frameIdx)+frameByteIdx*8);
            %if this is a normal sync word, use normal bits
            for bit_idx=0:7  %bytes are 8 bits long ;)                            
                if(dataStreamIn(SyncWordAllIndex(frameIdx)+frameByteIdx*8+bit_idx)=='0')               
                    byte = bitshift(byte,1); %This is a zero, just shift           
                else                    
                    byte = bitshift(byte,1); %This is a one, set the bit then shift               
                    byte = bitor(byte,1);              
                end                
            end        
        minorFrames(frameIdx,frameByteIdx+1)=byte;    
        end
    else %this minor frame is inverted
        for frameByteIdx=0:103 %minor frames are 103 bytes long
            byte=0;
            %Start of byte time
            frameTime(frameIdx,frameByteIdx+1)=bitTime(SyncWordAllIndex(frameIdx)+frameByteIdx*8);            
            for bit_idx=0:7  %bytes are 8 bits long ;)                            
                if(dataStreamIn(SyncWordAllIndex(frameIdx)+frameByteIdx*8+bit_idx)=='0')               
                    byte = bitshift(byte,1); %This is a zero, just shift
                    byte = bitor(byte,1);              
                else                    
                    byte = bitshift(byte,1); %This is a one, set the bit then shift                                  
                end                
            end        
        minorFrames(frameIdx,frameByteIdx+1)=byte;    
        end
    end
end
end