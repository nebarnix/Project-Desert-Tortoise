
oldest = 1;
framesFound = 0;
zero=0;
one=1;
byte=0;
syncWordLength = 19;
historyBufferCirc = zeros(syncWordLength,1)+48;
minorFrameShiftFlag = 0;
bitIdx = 1;
syncWord = '1110110111100010000';
%syncWord = '0001001000011101111';

%first one starts at 515
for idx = 1:numel(manchesterStreamBits) 
%for idx = (514):1400
%for idx = (514):1600

    %Enter this loop if we found a syncword last time around
    if(minorFrameShiftFlag == 1)
        
        if(manchesterStreamBits(idx)=='0')
            byte = bitshift(byte,1); %byte = byte << 1; %This is a zero, just shift
            byte = bitor(byte,zero); %byte = byte | zero;            
        else
            byte = bitshift(byte,1); % byte = byte << 1; %This is a one, set the bit then shift
            byte = bitor(byte,one); %byte = byte | one;
        end
        
        bitIdx = bitIdx + 1;
        
        if(bitIdx > 8)
            %minorFrame[frameByteIdx]=byte;
            fprintf('%.2X ',byte);
            byte = 0;
            bitIdx = 1;
            frameByteIdx = frameByteIdx + 1;
            
            if(frameByteIdx > 103)
                minorFrameShiftFlag = 0;
                fprintf('\n');
            end
        end
    end
    
    %overwrite oldest bit in cir buffer with newest bit
    historyBufferCirc(oldest+1) = manchesterStreamBits(idx);
    %fprintf('h[%d]=%c\n',oldest+1, historyBufferCirc(oldest+1));
    syncIndicator = 1;
    
    
    %Look for syncword
    for idx2 = 1:syncWordLength
        %compare syncword bytes to appropriate circular buffer bytes
        %fprintf('%d ',mod((oldest + idx2 + 1), syncWordLength)+1);
        if (syncWord(idx2) ~= historyBufferCirc(mod((oldest + idx2 + 1), syncWordLength)+1))
            syncIndicator = 0;
            break;
        end
    end
    %fprintf('\n');
    
    if(syncIndicator == 1 && minorFrameShiftFlag == 0)
        %fprintf('%.5f ',manchesterStreamBitsTime(idx));
        fprintf('%.2X ', hex2dec('ED'));
        fprintf('%.2X ', hex2dec('E2'));
        frameByteIdx = 2;
        minorFrameShiftFlag = 1;
        framesFound = framesFound+1;
        bitIdx = 4;
        byte = 0;
        zero = 0;
        one = 1;
    end
    
    
    %Look for Inverse Syncword
    syncIndicator = 1;
    for idx2 = 1:syncWordLength        
        %compare syncword bytes to appropriate circular buffer bytes
        %fprintf('sw[%d] == hst[%d] , %c == %c\n', idx2,mod((oldest + idx2) , syncWordLength),syncWord(idx2),historyBufferCirc(mod((oldest + idx2), syncWordLength)+1));
        %fprintf('%d ',mod((oldest + idx2 + 1), syncWordLength)+1);
        if syncWord(idx2) == historyBufferCirc(mod((oldest + idx2), syncWordLength)+1)
            syncIndicator = 0;
            %fprintf('%d==%d\n',historyBufferCirc(mod((oldest + idx2), syncWordLength)+1), syncWord(idx2));
            break;
        end
    end
    %fprintf('\n');
    
    if(syncIndicator == 1 && minorFrameShiftFlag == 0)
        
        %fprintf('%.5fi ',manchesterStreamBitsTime(idx));
        fprintf('%.2Xi ', hex2dec('ED'));
        fprintf('%.2X ', hex2dec('E2'));
        frameByteIdx = 2;
        minorFrameShiftFlag = 1;
        framesFound = framesFound + 1;
        bitIdx=4;
        byte=0;
        zero = 1;
        one = 0;
    end
    
    if minorFrameShiftFlag == 2
        for idx2 = 1:syncWordLength
            %compare syncword bytes to appropriate circular buffer bytes
            %fprintf('sw[%d] == hst[%d] , %c == %c\n', idx2,mod((oldest + idx2) , syncWordLength),syncWord(idx2),historyBufferCirc(mod((oldest + idx2), syncWordLength)+1));
            fprintf('%c',historyBufferCirc(mod((oldest + idx2 ), syncWordLength)+1));
        end
        fprintf('\n');
        
        for idx2 = 1:syncWordLength
            fprintf('%c',syncWord(idx2));
        end
        fprintf('\n');
        fprintf('\n');
    end
    
    %advance oldest bit pointer
    oldest = mod((oldest + 1), syncWordLength);
end

fprintf('%d Frames Found\n',framesFound);

%syncWordDetect(manchesterStreamBits(514:3344));

