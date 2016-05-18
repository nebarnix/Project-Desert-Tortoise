%static gain block makes the average value of the data set = 1
%Takes n points from the input
function [gainOut] = StaticGain(dataStreamIn,n,desiredGain)
%n=100000; %number of points to average
if nargin < 2 || n==0 
    n = numel(dataStreamIn);
end
avgSignalStrength = sum(abs(dataStreamIn(1:n)))/n;
gainOut = desiredGain/avgSignalStrength;
%dataStreamIn = gain.*dataStreamIn;
end