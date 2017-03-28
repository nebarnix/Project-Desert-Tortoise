for idx = 1:size(minorFrames,1)
    for idx2 = 1:size(minorFrames,2)
        fprintf('%0.2X ',minorFrames(idx,idx2));
    end
    fprintf('\n');
end