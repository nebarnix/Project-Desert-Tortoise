%Automatic Gain Control Block
function [dataStreamOut, gaini] = NormalizingAGC(dataStreamIn,AGC_loop_gain)
%Todo: implement a 'relock' mode for LARGE error values (either adjust gain
%outright or adjust loop gain

dataStreamOut = zeros(1,numel(dataStreamIn));
gaini = zeros(1,numel(dataStreamIn));

gain = 1; %Initial Gain Value
desired = 0.6366; %because a sin wave of amplitude 1 has this average absolute value
%AGC_loop_gain = ; 

progress = 0;
percentcomplete=0;
onePercent = numel(dataStreamOut) / 100;
msgCount = 1;
fprintf('Normalizing AGC:');
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
    
    dataStreamOut(idx) = gain * dataStreamIn(idx);
    error = desired - (gain * abs(dataStreamIn(idx)));
    
    %if(error > 1.0) %do something for really large errors
    %    error = -1/error;
    %end
    
    gain = gain + AGC_loop_gain * error ;
    gaini(idx) = gain;
end
for progress=1:msgCount-1      
    fprintf(char(8));        
end
fprintf('100%%\n');
end