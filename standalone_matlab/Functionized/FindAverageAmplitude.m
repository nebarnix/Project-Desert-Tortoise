function [ dataStreamOut ] = FindAverageAmplitude( dataStreamIn, alpha )
%FINDAVERAGEAMPLITUDE Find the smooth amplitude of a signal
%   returns a moving window average of the signal
%amplitude
dataStreamOut = zeros(1,numel(dataStreamIn));

average = 0;
progress = 0;
percentcomplete=0;
onePercent = numel(dataStreamOut) / 100;
msgCount = 1;
fprintf('Computing Raw Average:');
MSG = [num2str(percentcomplete) '%%'];
        msgCount = numel(MSG);
        fprintf(MSG);

for idx=1:numel(dataStreamIn)
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
    
    average = average * (1.0 - alpha) + alpha*abs(dataStreamIn(idx));
    dataStreamOut(idx) = average;
end

for progress=1:msgCount-1      
    fprintf(char(8));        
end
fprintf('100%%\n');

end

