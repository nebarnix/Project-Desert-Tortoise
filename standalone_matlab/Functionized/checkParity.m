%PARITY - Run Partiy Check on frames
function [goodFrames, parity] = checkParity(minorFrames)
%Word 103
%Bit 1: CPU B data transfer incomplete bit
%Bit 2: CPU A data transfer incomplete bit
%Bit 3: Even parity check in words 2 through 18
%Bit 4: Even parity check in words 19 thru 35
%Bit 5: Even parity check in words 36 thru 52
%Bit 6: Even parity check in words 53 thru 69
%Bit 7: Even parity check in words 70 thru 86
%Bit 8: Even parity check in words 87 thru bit 7 of word 103

%Count ones. If the number of one's is odd, modulus will be 1. Even
%parity bit will be set to change the number of ones to be even, so the
%parity bit should match the modulus of the one count exactly. 
%clear goodFrames parity;
parity = zeros(size(minorFrames,1),5);
goodFrames = 0;

for frame=1:size(minorFrames,1)
    parity(frame,1)=0;
    for word=3:19
        byte = minorFrames(frame,word);
           for shift=0:7
              parity(frame,1) = parity(frame,1)+bitand(bitshift(byte,-shift),1);   
           end       
    end
    
    for word=20:36           
        byte = minorFrames(frame,word);
           for shift=0:7
              parity(frame,2) = parity(frame,2)+bitand(bitshift(byte,-shift),1);   
           end       
    end
    
    for(word=37:53)           
        byte = minorFrames(frame,word);
           for(shift=0:7)
              parity(frame,3) = parity(frame,3)+bitand(bitshift(byte,-shift),1);   
           end       
    end
    
    for(word=54:70)           
        byte = minorFrames(frame,word);
           for(shift=0:7)
              parity(frame,4) = parity(frame,4)+bitand(bitshift(byte,-shift),1);   
           end       
    end
    
    for word=71:87
        byte = minorFrames(frame,word);
           for(shift=0:7)
              parity(frame,5) = parity(frame,5)+bitand(bitshift(byte,-shift),1);   
           end       
    end
        
    if(mod(parity(frame,1),2) == bitand(bitshift(minorFrames(frame,104),-5),1)) %check if divisible by 2 (even)       
        parity(frame,1) = 0; %words might be good or might have an even number of bit errors         
    else        
        parity(frame,1) = 1; %Words contain at least one error
    end    
        
    if(mod(parity(frame,2),2) == bitand(bitshift(minorFrames(frame,104),-4),1)) %check if divisible by 2 (even)       
        parity(frame,2) = 0; %words might be good or might have an even number of bit errors         
    else        
        parity(frame,2) = 1; %Words contain at least one error
    end
    
    if(mod(parity(frame,3),2) == bitand(bitshift(minorFrames(frame,104),-3),1)) %check if divisible by 2 (even)       
        parity(frame,3) = 0; %words might be good or might have an even number of bit errors         
    else        
        parity(frame,3) = 1; %Words contain at least one error
    end
    
    if(mod(parity(frame,4),2) == bitand(bitshift(minorFrames(frame,104),-2),1)) %check if divisible by 2 (even)       
        parity(frame,4) = 0; %words might be good or might have an even number of bit errors         
    else        
        parity(frame,4) = 1; %Words contain at least one error
    end
    
    if(mod(parity(frame,5),2) == bitand(bitshift(minorFrames(frame,104),-1),1)) %check if divisible by 2 (even)       
        parity(frame,5) = 0; %words might be good or might have an even number of bit errors         
    else        
        parity(frame,5) = 1; %Words contain at least one error
    end
    
    if sum(parity(frame,:)) == 0
        goodFrames = goodFrames + 1;
    end
end
fprintf(['\n' num2str(goodFrames) ' out of ' num2str(frame) ' Error Free Frames\n\n']);
fprintf([num2str(numel(parity(parity == 0))) ' Good Chunks and ' num2str(numel(parity(parity == 1))) ' Bad Chunks\n\n']);
