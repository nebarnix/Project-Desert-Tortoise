%BITSTREAM - Look for syncword and its inverse (in case of phase reversal)
function [SyncWordIndex, SyncWordInvIndex] = syncWordDetect(dataStreamIn)

SyncWord = '0001011110000';
SyncWordInverse = '1110100001111';

%S1 = '11101101111000100000 %1101'; %NOAA18 ID last 4 1101
%S2 = '00010010000111011111 0010'; %NOAA18 ID

SyncWordIndex = strfind(dataStreamIn, SyncWord)-3;
SyncWordInvIndex = strfind(dataStreamIn, SyncWordInverse)-3;
fprintf([ '\n' num2str(numel(SyncWordInvIndex)+numel(SyncWordIndex)) ' detected\n']);

%plot(k1(2:end),diff(k1),'o',k2(2:end),diff(k2),'x',errx,erry*10,'.');
%plot(k1,1:length(k1),'o',k2,1:length(k2),'x',errx,erry*10,'.');

%If you want plots, uncomment
%figure(3);
%subplot(2,1,1);
%plot(bitTime(SyncWordIndex),SyncWordIndex,'o',bitTime(SyncWordInvIndex),SyncWordInvIndex,'x',bitTime((errx(1:end-1))),erry(1:(end-1))*10,'.');
%subplot(2,1,2);
%plot(bitTime(SyncWordIndex(2:end)),diff(SyncWordIndex),'o',bitTime(SyncWordInvIndex(2:end)),diff(SyncWordInvIndex),'x',bitTime(errx(1:(end-1))),erry(1:(end-1))*10,'.');

