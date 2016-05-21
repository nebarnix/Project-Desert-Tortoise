function [dataStreamOut] = Squelch(dataStreamIn, squelchSigIn, squelchThreshold)
dataStreamOut = zeros(1,numel(dataStreamIn));

for idx=1:numel(dataStreamIn)
    if squelchSigIn(idx) > squelchThreshold
        dataStreamOut(idx) = dataStreamIn(idx);        
    end
end
