%BITSTREAM - Look for syncword and its inverse (in case of phase reversal)
function [SyncWordIndex, SyncWordInvIndex] = syncWordDetect(dataStreamIn)

SyncWord = '1110110111100010000';
SyncWordInverse = '0001001000011101111';

%S1 = '11101101111000100000 %1101'; %NOAA18 ID last 4 1101
%S2 = '00010010000111011111 0010'; %NOAA18 ID

SyncWordIndex = strfind(dataStreamIn, SyncWord);
SyncWordInvIndex = strfind(dataStreamIn, SyncWordInverse);
fprintf([ '\n' num2str(numel(SyncWordInvIndex)+numel(SyncWordIndex)) ' detected\n' num2str(sum(mod(diff(SyncWordIndex),832)==0) + sum(mod(diff(SyncWordInvIndex),832)==0)) ' match length\n\n']);

%plot(k1(2:end),diff(k1),'o',k2(2:end),diff(k2),'x',errx,erry*10,'.');
%plot(k1,1:length(k1),'o',k2,1:length(k2),'x',errx,erry*10,'.');

%If you want plots, uncomment
%figure(3);
%subplot(2,1,1);
%plot(bitTime(SyncWordIndex),SyncWordIndex,'o',bitTime(SyncWordInvIndex),SyncWordInvIndex,'x',bitTime((errx(1:end-1))),erry(1:(end-1))*10,'.');
%subplot(2,1,2);
%plot(bitTime(SyncWordIndex(2:end)),diff(SyncWordIndex),'o',bitTime(SyncWordInvIndex(2:end)),diff(SyncWordInvIndex),'x',bitTime(errx(1:(end-1))),erry(1:(end-1))*10,'.');

