%% *LOAD - load raw IQ data and prepare for decoding
%
clear all;
%clf;
%hfile = 'POES_56k250.raw';
%hfile = 'pll_nocarrier_polyphased_equalized_loproto2.raw';
%hfile = 'C:\Users\nebarnix\Documents\vmshare\POES\6-NOAA15\6_bits_dcblock.raw'; 
%hfile = 'C:\Users\nebarnix\Documents\vmshare\POES\38-N18 SDRplay\bits_grc.raw';
hfile = 'C:\Users\nebarnix\Documents\vmshare\POES\33-N18 hackrf\33_bits_50k.raw';
%hfile = 'e:\N15_Aug0415_bits.raw';
%hfile = 'e:\N18_Aug0415_bits.raw';
%hfile = '06-28-39_137350kHz_bits.raw'; 
fid = fopen(hfile,'rb');
y = fread(fid,'float32');

y = (y(1:2:end) + 1i*y(2:2:end))';
%y = (1 + 1i*y(1:end))';
%y = y(2:2:end) + 1i*y(1:2:end);
Nfft = 1500000;
%GNUradio should generate 1 sample per symbol with a period of
%- (1 /(8320*2)) Seconds per sample
Ts = 1/(8320*2);
RawTime=0:Ts:Ts*(numel(y)-1);
fclose(fid);
%% *Import data from MM script
clear y;
y=imag(dataStreamBits)+1i*real(dataStreamBits); %quick swap of real and imag data
Ts = 1/(8320*2);
%RawTime=0:Ts:Ts*(numel(y)-1);
RawTime=dataStreamBitsTime;

%% CONSTELLATION - plot bit strengths
pointScale = 51;
length(y(1:pointScale:end));
scatter(RawTime(1:pointScale:end),imag(y(1:pointScale:end)),2,'.');
axis([0 max(RawTime) -2 2]);

%% CONSTELLATION - Show 10000 points again and again as run
%the constellation should be two points VERTICALLY ALIGNED else swap I/Q above
figure(1); %164 seconds set 6, best
colormap copper;
%Nfft = 100000+Nfft
x = y(Nfft:10000+Nfft);
c = linspace(100,1,length(x));
h=scatter(real(x),imag(x),2,c);
xlabel('real');
ylabel('imag');
axis([-2 2 -2 2]);
axis square;
%text(-1.8,1.8,['T=' num2str(RawTime(Nfft),4) '-' num2str(RawTime(10000+Nfft),4)]);
title(['Signal Constellation at time interval T=' num2str(RawTime(Nfft),4) '-' num2str(RawTime(10000+Nfft),4) 's']);

%% CONSTELLATION - Animate
figure(1);
colormap copper;
scatter(1,1);

axis([-2 2 -2 2]);
numPoints=1000;
c = linspace(1,10,numPoints+1);

for Nfft=1:numPoints:numel(y)-numPoints
    x = y(Nfft:numPoints+Nfft);
    scatter(real(x),imag(x),[],c);
    axis([-2 2 -2 2]);
    text(-1.8,1.8,num2str(RawTime(Nfft)));
    %grid on;
    %Nfft
    drawnow
end

%% *CONSTELLATION - convert to bits from raw manchester bits
%resyncThreshold = 0.75; %set 6 - 3599 out of 4066 Error Free Frames
%resyncThreshold = 0.80; %set 6 - 3602 out of 4067 Error Free Frames
%resyncThreshold = 0.85; %set 6 - 3602 out of 4071 Error Free Frames *sweet spot for set 6
%resyncThreshold = 0.9; %set 6 -- 3602 out of 4071 Error Free Frames
%resyncThreshold = 1.0; %set 6 -- 3601 out of 4070 Error Free Frames
resyncThreshold = 1;
%resyncThreshold = 0.50; %set 20 - 13872 Good Chunks and 6293 Bad Chunks 1110 out of 4033 Error Free Frames
%******THIS ONE %resyncThreshold = 0.60; %set 20 - 14154 Good Chunks and 6146 Bad Chunks 1142 out of *4060 Error Free Frames
%resyncThreshold = 0.65; %set 20 - 14123 Good Chunks and 6042 Bad Chunks 1152 out of 4033 Error Free Frames
%resyncThreshold = 0.70; %set 20 - 14161 Good Chunks and 5909 Bad Chunks 1166 out of 4014 Error Free Frames
%resyncThreshold = 0.75; %set 20 - *14209 Good Chunks and 5786 Bad Chunks 1164 out of 3999 Error Free Frames
%resyncThreshold = 0.80; %set 20 - 14135 Good Chunks and 5685 Bad Chunks 1165 out of 3964 Error Free Frames
%resyncThreshold = 0.85; %set 20 - 14063 Good Chunks and 5602 Bad Chunks *1168 out of 3933 Error Free Frames
%resyncThreshold = 0.90; %set 20 - 14021 Good Chunks and 5614 Bad Chunks 1162 out of 3927 Error Free Frames

%Idea - this threshold could be made dynamic and could be a percentage of
%a moving average over the last X points. This would let the algorithm
%follow gain->amplitude variations caused by interference from nearby signals. 

idx2=1;
clockmod = 0;
idxerr=1;
%idxerr2=1;
clear errx erry bitstreamManchester BitTime;
bitstreamManchester(1) = uint8(0);
for idx = 2:numel(y)-1    
    %If not a bit boundary, see if it should be and we're out of sync
    %But only resync on strong bits
    if(mod(idx,2) ~= clockmod)    
        if(sign(imag(y(idx-1))) == sign(imag(y(idx))))
            errx(idxerr)=idx2;
            erry(idxerr)=imag(y(idx));
            idxerr=idxerr+1;   
            if(abs(imag(y(idx-1))) > resyncThreshold && abs(imag(y(idx))) > resyncThreshold)                
                clockmod = mod(idx,2); %only resync if we have confidence in the bit decisions
            end
            
        end        
    end
    
    %check for bit boundary, and make decision using the strongest of the
    %two bits. 
    if(mod(idx,2) == clockmod)
        if(abs(imag(y(idx))) > abs(imag(y(idx+1)))) %use the strongest symbol to determine bit
            if(imag(y(idx)) > 0)                
                bitstreamManchester(idx2) = '0';
            else
                bitstreamManchester(idx2) = '1';
            end
        else
            if(imag(y(idx+1)) > 0)                
                bitstreamManchester(idx2) = '1';
            else
                bitstreamManchester(idx2) = '0';
            end        
        end
        
        BitTime(idx2)=RawTime(idx);    
        idx2 = idx2+1;
        %bitstream_manchester(idx2) = bitstream(idx+1);
                               
        
        %if(bitstream(idx-1) == bitstream(idx))
        %    errx(idxerr)=idx2;
        %    erry(idxerr)=real(y(idx));
        %    idxerr=idxerr+1;
        %end            
    end        
    
    %look for two same bits in a row
    %this can only happen on bit boundary
    %if we find this, resync to this boundary
    %if(bitstream(idx-1) == bitstream(idx))
    %    if(bitstream(idx) == 0)
    %        clockmod = mod(idx,2);
    %    end
    %end
end

%%oldcode below
%%plot autocorrelation to see if anything is here
%figure(2);
%subplot(2,1,1);
%plot(xcorr(bitstream))
%subplot(2,1,2);
%plot(xcorr(bitstream_manchester))

%%test autocorrelation of random number generator as reference
%testarray = rand(1,100000);
%for idx = 1:numel(testarray)
%    if(idx>0.5)
%        testarray2(idx) = 1;
%    else
%        testarray2(idx) = 0;
%    end
%end
%plot(xcorr(testarray2))

%% *BITSTREAM - Look for syncword and its inverse (in case of phase reversal)
SyncWord = '1110110111100010000'; %0100'; %NOAA15 ID last 4 0100
SyncWordInverse = '0001001000011101111'; %1011'; %NOAA15 ID

%S1 = '11101101111000100000 %1101'; %NOAA18 ID last 4 1101
%S2 = '00010010000111011111 0010'; %NOAA18 ID

SyncWordIndex = strfind(bitstreamManchester, SyncWord);
SyncWordInvIndex = strfind(bitstreamManchester, SyncWordInverse);
fprintf([ '\n' num2str(numel(SyncWordInvIndex)+numel(SyncWordIndex)) ' detected\n' num2str(sum(mod(diff(SyncWordIndex),832)==0) + sum(mod(diff(SyncWordInvIndex),832)==0)) ' match length\n\n' num2str(idxerr) ' errors\n\n']);

%plot(k1(2:end),diff(k1),'o',k2(2:end),diff(k2),'x',errx,erry*10,'.');
%plot(k1,1:length(k1),'o',k2,1:length(k2),'x',errx,erry*10,'.');

%If you want plots, uncomment
%figure(3);
%subplot(2,1,1);
%plot(BitTime(SyncWordIndex),SyncWordIndex,'o',BitTime(SyncWordInvIndex),SyncWordInvIndex,'x',BitTime((errx(1:end-1))),erry(1:(end-1))*10,'.');
%subplot(2,1,2);
%plot(BitTime(SyncWordIndex(2:end)),diff(SyncWordIndex),'o',BitTime(SyncWordInvIndex(2:end)),diff(SyncWordInvIndex),'x',BitTime(errx(1:(end-1))),erry(1:(end-1))*10,'.');

%% *BITSTREAM -  convert ascii binary stream to actual binary at matched syncword locations
%Includes correct and inverted syncwords for periods of constellation reversal. 

clear minorFrames frameTime;

SyncWordAllIndex = sort(cat(2,SyncWordIndex,SyncWordInvIndex));

for frameIdx=1:numel(SyncWordAllIndex)-1
    %See if the frame is normal or inverted bits
    if isempty(find(SyncWordInvIndex == SyncWordAllIndex(frameIdx),1))
        for frameByteIdx=0:103 %minor frames are 103 bytes long
            byte=0;
            %Start of byte time
            frameTime(frameIdx,frameByteIdx+1)=BitTime(SyncWordAllIndex(frameIdx)+frameByteIdx*8);
            %if this is a normal sync word, use normal bits
            for bit_idx=0:7  %bytes are 8 bits long ;)                            
                if(bitstreamManchester(SyncWordAllIndex(frameIdx)+frameByteIdx*8+bit_idx)=='0')               
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
            frameTime(frameIdx,frameByteIdx+1)=BitTime(SyncWordAllIndex(frameIdx)+frameByteIdx*8);            
            for bit_idx=0:7  %bytes are 8 bits long ;)                            
                if(bitstreamManchester(SyncWordAllIndex(frameIdx)+frameByteIdx*8+bit_idx)=='0')               
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
%figure(4);
%plot(bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6))); %plot minor frame counter, which is annoyingly 9 bits long
%plot(frameTime(:,5),bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6)),'.');

%% *BITSTREAM - Print Major Frame UTC times and spacecraft ID
clear minorFrameID dayNum spacecraftID;
minorFrameID = bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6)); %pull out minor frame counter, which is annoyingly 9 bits long

%If you want plots, uncomment below
%figure(4);
%plot(frameTime(:,5),minorFrameID,'.'); %plot minor frame counter

spaceCraft(1) = 0;
dayNum(1) = 0;

idx = 1;

for frame=1:numel(minorFrameID)        
    spaceCraft(frame) = minorFrames(frame,3);    
    
    if(minorFrameID(frame) == 0)
        dayNum(idx)=bitshift(minorFrames(frame,8+1),1)+bitshift(bitor(minorFrames(frame,9+1),128),-7); %should be 241 for loproto2 recording        
        
        %check if the counter is more than the number of ms per day (bad data)
        if((bitshift(bitand(minorFrames(frame,9+1),bin2dec('111')),24) + bitshift(minorFrames(frame,10+1),16) + bitshift(minorFrames(frame,11+1),8) + minorFrames(frame,12+1)) < 86400000)
            dayMSeconds(idx) = (bitshift(bitand(minorFrames(frame,9+1),bin2dec('111')),24) + bitshift(minorFrames(frame,10+1),16) + bitshift(minorFrames(frame,11+1),8) + minorFrames(frame,12+1));
            fprintf([num2str(frameTime(frame,1)) ' Local Seconds is ' num2str(dayMSeconds(idx)/1000.0) ' Spacecraft Day Seconds' ]);
            
            hour = floor(dayMSeconds(idx)/(1000*60*60));
            minute = floor((dayMSeconds(idx)/(1000*60*60) - hour)*60);
            seconds = (((dayMSeconds(idx)/(1000*60*60) - hour)*60) - minute) *60;
            fprintf([' which is ' num2str(hour) ':' num2str(minute) ':' num2str(seconds)]);
            
            if idx > 1 && dayMSeconds(idx) <= dayMSeconds(idx-1)
                fprintf(' ...but this might be an error');
            else
                hour = floor((dayMSeconds(idx)-frameTime(frame,1)*1000)/(1000*60*60));
                minute = floor(((dayMSeconds(idx)-frameTime(frame,1)*1000)/(1000*60*60) - hour)*60);
                seconds = ((((dayMSeconds(idx)-frameTime(frame,1)*1000)/(1000*60*60) - hour)*60) - minute) *60;
                fprintf(['    t(0) => ' num2str(hour) ':' num2str(minute) ':' num2str(seconds)]);
            end
            fprintf('\n');
        else
            dayMSeconds(idx) = -1;
        end        
        idx = idx + 1;
        %fprintf([num2str(frameTime(frame,1)) ':' num2str(dec2bin(minorFrames(frame,9+1),8)) ' ' num2str(dec2bin(minorFrames(frame,10+1),8)) ' ' num2str(dec2bin(minorFrames(frame,11+1),8)) ' ' num2str(dec2bin(minorFrames(frame,12+1),8)) '\n']);
        %fprintf([num2str(day) ' ' num2str(dec2hex(minorFrames(frame,8+1))) ' ' num2str(dec2hex(minorFrames(frame,9+1))) '\n\n']);        
    end
end

if mode(spaceCraft) == 8
   fprintf(['Spacecraft: ' num2str(mode(spaceCraft)) '=>NOAA-15\n']);    
   fprintf(['Julean Day: ' num2str(mode(dayNum)) ' \n']);
elseif mode(spaceCraft) == 13
    fprintf(['Spacecraft: ' num2str(mode(spaceCraft)) '=>NOAA-18\n']);    
    fprintf(['Julean Day: ' num2str(mode(dayNum)) ' \n']);
elseif mode(spaceCraft) == 15
    fprintf(['Spacecraft: ' num2str(mode(spaceCraft)) '=>NOAA-19\n']);    
    fprintf(['Julean Day: ' num2str(mode(dayNum)) ' \n']);
else
    fprintf(['Spacecraft: ' num2str(mode(spaceCraft)) '=>A UFO!\n']);    
end


%% *PARITY - Run Partiy Check on frames
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
clear goodFrames parity;
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
fprintf([num2str(goodFrames) ' out of ' num2str(frame) ' Error Free Frames\n\n']);
fprintf([num2str(numel(parity(parity == 0))) ' Good Chunks and ' num2str(numel(parity(parity == 1))) ' Bad Chunks\n\n']);

%% PARITY - Plot results
figure(5);
clf;
%plot(1:size(parity,1),parity(:,1),'.-');
%plot(bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6)));
plot(frameTime(:,5),bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6)),'.');
%Px=1:size(parity,1);
Px=frameTime(:,1)';
Py=parity(:,1)'*512;
Pz = zeros(1,size(Px,2));
S = surface([Px;Px],[Py;Py],[Pz;Pz], 'facecol','no','edgecol','interp','linew',1.5,'edgealpha',.125,...
            'edgecolor',0*[1 1 1],...
            'facecolor',0*[1 1 1]);

%Px=1:size(parity,1);
%Px=frameTime(:,1)';
Py=parity(:,2)'*512;
Pz = zeros(1,size(Px,2));
S = surface([Px;Px],[Py;Py],[Pz;Pz], 'facecol','no','edgecol','interp','linew',1.5,'edgealpha',.125,...
            'edgecolor',0.5*[.15 1 1],...
            'facecolor',0.5*[.15 1 1]);

%Px=1:size(parity,1);
Py=parity(:,3)'*512;
Pz = zeros(1,size(Px,2));
S = surface([Px;Px],[Py;Py],[Pz;Pz], 'facecol','no','edgecol','interp','linew',1.5,'edgealpha',.125,...
            'edgecolor',0.5*[1 .15 1],...
            'facecolor',0.5*[1 .15 1]);

 %Px=1:size(parity,1);
 Py=parity(:,4)'*512;
 Pz = zeros(1,size(Px,2));
 S = surface([Px;Px],[Py;Py],[Pz;Pz], 'facecol','no','edgecol','interp','linew',1.5,'edgealpha',.125,...
             'edgecolor',0.5*[1 1 .15],...
             'facecolor',0.5*[1 1 .15]);
 
%Px=1:size(parity,1);
Py=parity(:,5)'*512;
Pz = zeros(1,size(Px,2));
S = surface([Px;Px],[Py;Py],[Pz;Pz], 'facecol','no','edgecol','interp','linew',1.5,'edgealpha',.125,...
            'edgecolor',0.5*[.15 1 .15],...
            'facecolor',0.5*[.15 1 .15]);

axis([0 max(frameTime(:,1)) -0.1 515]);

%% Telemetry
%STX1 Power, TIP 16-sec analog subcom-1, 11, 48,128,208,288 - high-gain HRPT antenna
%STX2 Power, TIP 16-sec analog subcom-1, 11, 50,130,210,290 - OMNI HRPT antenna
%STX3 Power, TIP 16-sec analog subcom-1, 11, 40,120,200,280 - 
%SARR-A Power, TIP 16-sec analog subcom-2, 14, 114,274
%SARR-B Power, TIP 16-sec analog subcom-2, 14, 2,162
clear AnalogA1TLM AnalogA2TLM AA1SARRAPwr AA1SARRBPwr AA1STX1Pwr AA1STX2Pwr AA1STX3Pwr;
clear AA1SARRAPwrT AA1SARRBPwrT AA1STX1PwrT AA1STX2PwrT AA1STX3PwrT;
clear titleString;

AA1SARRAPwr(1)=0;
AA1SARRBPwr(1)=0;
AA1STX1Pwr(1)=0;
AA1STX2Pwr(1)=0;
AA1STX3Pwr(1)=0;
AA1SARRAPwrT(1)=0;
AA1SARRBPwrT(1)=0;
AA1STX1PwrT(1)=0;
AA1STX2PwrT(1)=0;
AA1STX3PwrT(1)=0;

AnalogA1TLM(:,1) = minorFrames(:,11+1);
AnalogA1TLM(:,2) = frameTime(:,11+1);
AnalogA1TLM(:,3) = minorFrameID(:);

AnalogA2TLM(:,1) = minorFrames(:,14+1);
AnalogA2TLM(:,2) = frameTime(:,14+1);
AnalogA2TLM(:,3) = minorFrameID(:);

for idx=1:size(AnalogA1TLM,1)
    if(AnalogA1TLM(idx,3) == 48 || AnalogA1TLM(idx,3) == 128 || AnalogA1TLM(idx,3) == 208 || AnalogA1TLM(idx,3) == 288)
        AA1STX1Pwr(end+1) = AnalogA1TLM(idx,1);    
        AA1STX1PwrT(end+1) = AnalogA1TLM(idx,2);
    elseif(AnalogA1TLM(idx,3) == 50 || AnalogA1TLM(idx,3) == 130 || AnalogA1TLM(idx,3) == 210 || AnalogA1TLM(idx,3) == 290)
        AA1STX2Pwr(end+1) = AnalogA1TLM(idx,1);
        AA1STX2PwrT(end+1) = AnalogA1TLM(idx,2);
    elseif(AnalogA1TLM(idx,3) == 40 || AnalogA1TLM(idx,3) == 120 || AnalogA1TLM(idx,3) == 200 || AnalogA1TLM(idx,3) == 280)
        AA1STX3Pwr(end+1) = AnalogA1TLM(idx,1);
        AA1STX3PwrT(end+1) = AnalogA1TLM(idx,2);
    end
    
    if(AnalogA2TLM(idx,3) == 114 || AnalogA2TLM(idx,3) == 274)
        AA1SARRAPwr(end+1) = AnalogA2TLM(idx,1);
        AA1SARRAPwrT(end+1) = AnalogA2TLM(idx,2);
    elseif(AnalogA1TLM(idx,3) == 2 || AnalogA2TLM(idx,3) == 162)
        AA1SARRBPwr(end+1) = AnalogA2TLM(idx,1);
        AA1SARRBPwrT(end+1) = AnalogA2TLM(idx,2);
    end
end
plot(AA1STX1PwrT,AA1STX1Pwr,AA1STX2PwrT,AA1STX2Pwr,AA1STX3PwrT,AA1STX3Pwr,AA1SARRAPwrT,AA1SARRAPwr,AA1SARRBPwrT,AA1SARRBPwr);
titleString = sprintf('STX1=%d, STX2=%d, STX3=%d, SARR-A=%d, SARR-B=%d ',mode(AA1STX1Pwr),mode(AA1STX2Pwr),mode(AA1STX3Pwr),mode(AA1SARRAPwr),mode(AA1SARRBPwr));;
title(titleString);

%% HIRS - prepare HIRS data structures
clear HIRSdata HIRSTime;
%Pull out HIRS3 data, which is actually an embedded stream of 13-bit words
%16,17,
%22,23,
%26,27,
%30,31,
%34,35,
%38,39,
%42,43,
%54,55,
%58,59,
%62,63,
%66,67,
%70,71,
%74,75,
%78,79,
%82,83,
%84,85,
%88,89,
%92,93

HIRSdata(:,1) = (minorFrames(:,16+1));
HIRSdata(:,2) = (minorFrames(:,17+1));
HIRSTime(:,1) = (frameTime(:,16+1));
HIRSTime(:,2) = (frameTime(:,17+1));

HIRSdata(:,3) = (minorFrames(:,22+1));
HIRSdata(:,4) = (minorFrames(:,23+1));
HIRSTime(:,3) = (frameTime(:,22+1));
HIRSTime(:,4) = (frameTime(:,23+1));

HIRSdata(:,5) = (minorFrames(:,26+1));
HIRSdata(:,6) = (minorFrames(:,27+1));
HIRSTime(:,5) = (frameTime(:,26+1));
HIRSTime(:,6) = (frameTime(:,27+1));

HIRSdata(:,7) = (minorFrames(:,30+1));
HIRSdata(:,8) = (minorFrames(:,31+1));
HIRSTime(:,7) = (frameTime(:,30+1));
HIRSTime(:,8) = (frameTime(:,31+1));

HIRSdata(:,9) = (minorFrames(:,34+1));
HIRSdata(:,10) = (minorFrames(:,35+1));
HIRSTime(:,9) = (frameTime(:,34+1));
HIRSTime(:,10) = (frameTime(:,35+1));

HIRSdata(:,11) = (minorFrames(:,38+1));
HIRSdata(:,12) = (minorFrames(:,39+1));
HIRSTime(:,11) = (frameTime(:,38+1));
HIRSTime(:,12) = (frameTime(:,39+1));

HIRSdata(:,13) = (minorFrames(:,42+1));
HIRSdata(:,14) = (minorFrames(:,43+1));
HIRSTime(:,13) = (frameTime(:,42+1));
HIRSTime(:,14) = (frameTime(:,43+1));

HIRSdata(:,15) = (minorFrames(:,54+1));
HIRSdata(:,16) = (minorFrames(:,55+1));
HIRSTime(:,15) = (frameTime(:,54+1));
HIRSTime(:,16) = (frameTime(:,55+1));

HIRSdata(:,17) = (minorFrames(:,58+1));
HIRSdata(:,18) = (minorFrames(:,59+1));
HIRSTime(:,17) = (frameTime(:,58+1));
HIRSTime(:,18) = (frameTime(:,59+1));

HIRSdata(:,19) = (minorFrames(:,62+1));
HIRSdata(:,20) = (minorFrames(:,63+1));
HIRSTime(:,19) = (frameTime(:,62+1));
HIRSTime(:,20) = (frameTime(:,63+1));

HIRSdata(:,21) = (minorFrames(:,66+1));
HIRSdata(:,22) = (minorFrames(:,67+1));
HIRSTime(:,21) = (frameTime(:,66+1));
HIRSTime(:,22) = (frameTime(:,67+1));

HIRSdata(:,23) = (minorFrames(:,70+1));
HIRSdata(:,24) = (minorFrames(:,71+1));
HIRSTime(:,23) = (frameTime(:,70+1));
HIRSTime(:,24) = (frameTime(:,71+1));

HIRSdata(:,25) = (minorFrames(:,74+1));
HIRSdata(:,26) = (minorFrames(:,75+1));
HIRSTime(:,25) = (frameTime(:,74+1));
HIRSTime(:,26) = (frameTime(:,75+1));

HIRSdata(:,27) = (minorFrames(:,78+1));
HIRSdata(:,28) = (minorFrames(:,79+1));
HIRSTime(:,27) = (frameTime(:,78+1));
HIRSTime(:,28) = (frameTime(:,79+1));

HIRSdata(:,29) = (minorFrames(:,82+1));
HIRSdata(:,30) = (minorFrames(:,83+1));
HIRSTime(:,29) = (frameTime(:,82+1));
HIRSTime(:,30) = (frameTime(:,83+1));

HIRSdata(:,31) = (minorFrames(:,84+1));
HIRSdata(:,32) = (minorFrames(:,85+1));
HIRSTime(:,31) = (frameTime(:,84+1));
HIRSTime(:,32) = (frameTime(:,85+1));

HIRSdata(:,33) = (minorFrames(:,88+1));
HIRSdata(:,34) = (minorFrames(:,89+1));
HIRSTime(:,33) = (frameTime(:,88+1));
HIRSTime(:,34) = (frameTime(:,89+1));

HIRSdata(:,35) = (minorFrames(:,92+1));
HIRSdata(:,36) = (minorFrames(:,93+1));
HIRSTime(:,35) = (frameTime(:,92+1));
HIRSTime(:,36) = (frameTime(:,93+1));

%% HIRS - Pull HIRS data into array
%288 bits mean 36 bytes or 22 words with 2 left over status bits (23 words)
%Easiest way I can think of in matlab is to convert to ASCII string of 1's
%and 0's like we did for the manchester decoding but there is probably (definately) a
%more efficient way
clear HIRSstream;
for frame=1:size(HIRSdata,1) %number of frames
    
    Intermediate = cellstr(dec2bin(HIRSdata(frame,:),8));
    HIRSstream(frame,:) = [Intermediate{:}];
    %HIRSstreamTime(frame,:) = HIRStime(frame,:)
    %HIRSstream(frame,:) = sprintf(dec2bin(HIRSdata(frame,:),8));
end

%% HIRS - Convert the strings of ascii bits to binary bits
%THIS DOES NOT GRAB THE LAST TWO BITS OF EACH STREAM YET (which
%are only status bits anyway ;)

%Words 3 through 22 contain a sign bit at the MSB which just means we can
%subtract the MSB value

%This is a good place to also calculate parity of the frame and compare to
%bit 288 - odd parity (this bit in inserted to maintain oddness!)

clear HIRSwords HIRSstreamParity;
HIRSstreamParity(length(HIRSdata)) = 0;
HIRSwords(1,1) = int16(0);
for frame=1:length(HIRSdata) %loop through each frame
   for idx=0:21 %loop through each minor frame's word
      byte = uint16(0);
      for bit_idx=0:12  %Loop through each words's bits     
           %idx*13+bit_idx+1
           if(HIRSstream(frame,idx*13+bit_idx+1)=='0')
               byte = bitshift(byte,1); %This is a zero, just shift
           else
               HIRSstreamParity(frame) = HIRSstreamParity(frame)+1; 
               %byte = bitor(byte,1);
               byte = bitshift(byte,1); %This is a one, set the bit then shift (first time sets leading zero, who cares)
               byte = bitor(byte,1);
           end
      end
      
      %Can't we do the same oepration simply by subtracting the MSB?
%       if(idx > 1)
%           HIRSwords(frame,idx+1) = int16(byte)-4096; %12 1's = 4095
%           %int16(byte)-4096
%       else          
%           HIRSwords(frame,idx+1) = uint16(byte);                                      
%       end
      
      if idx > 1 && bitshift(byte,-12) == 0 %if this is word 3-22, check the sign bit, then mask it out
          HIRSwords(frame,idx+1) = -int16(bitand(byte,4095)); %12 1's = 4095          
          %HIRSwords(frame,idx+1) = -4096;
          %fprintf([num2str(HIRSwords(frame,idx+1)) '\n']);
      elseif idx > 1 && bitshift(byte,-12) == 1  % positive, mask out sign bit
          HIRSwords(frame,idx+1) = int16(bitand(byte,4095)); %12 1's = 4095
      else          
           HIRSwords(frame,idx+1) = byte;                                      
      end      
   end
   
   %add in the data valid bit to the parity and calculate parity
   HIRSstreamParity(frame) = HIRSstreamParity(frame) + str2double(HIRSstream(frame,287));
   if(mod(HIRSstreamParity(frame),2) == 0 && HIRSstream(frame,288) == '1')
       HIRSstreamParity(frame) = 0; %good data, parity bit added to make even data odd
   else
       HIRSstreamParity(frame) = 1; %bad data
   end
end

%% HIRS - Plot raw data
%Plot scan Encoer position (increments of 1.8 degrees) and cal level
%indicator
ax1=subplot(3,1,1);
plot(HIRSTime(:,1),bitshift(uint16(HIRSwords(:,1)),-5),'x',HIRSTime(:,1),bitand(uint16(HIRSwords(:,1)),2^5-1),'o');
%plot(HIRSTime(:,1),uint16(HIRSwords(:,1)),'x');
title('Mirror Scan Position');

%Element number and period monitor (excpet I'm not sure what that is)
%0-55 are earth data
ax2=subplot(3,1,2);
plot(HIRSTime(:,2),bitand(bitshift(uint16(HIRSwords(:,2)),-1),63),'x',HIRSTime(:,2),bitshift(uint16(HIRSwords(:,2)),-7),'o'); %63 = bin2dec('111111')
title('element number');

ax3=subplot(3,1,3);
%plot(HIRSTime(:,8),HIRSwords(:,8),'x'); %63 = bin2dec('111111')
plot(HIRSTime(:,2),HIRSwords(:,3:end-1),'-x')
title('element X data');

linkaxes([ax1, ax2, ax3],'x');

%% HIRS - Separate Data by Filter and Plot Images
%There are 63 'elements'. 
%0-55 are earth view and each of the 20 words is an IR band
%56-63 are various parameter sets and the definition of the 20 words is
%specific to the set at hand.
clear testimage xdata ydata cdata;
close all;
maxElement = 64;

testimage(1,1) = uint16(0);

idx = 1;
idx2 = 1;

for frame=1:size(HIRSwords,1) %loop through each frame
   %HIRSElements() = bitand(bitshift(HIRSwords(:,2),-1),63);   
   Element = bitand(bitshift(uint16(HIRSwords(frame,2)),-1),63);   
   
   if sum(parity(frame,:)) == 0
      for(idx2 = 1:20) %associate all 20 words with a specific element
         cdataTrusted(idx,Element+1,idx2) = HIRSwords(frame,idx2+2)+4096;
      end    
   end
   %cdata(idx,Element+1) = HIRSwords(frame,3)+4096;
   
   for(idx2 = 1:20) %associate all 20 words with a specific element
      cdata(idx,Element+1,idx2) = HIRSwords(frame,idx2+2)+4096;
   end
   
%    cdata1(idx,Element+1) = HIRSwords(frame,3)+4096;
%    cdata2(idx,Element+1) = HIRSwords(frame,4)+4096;
%    cdata3(idx,Element+1) = HIRSwords(frame,5)+4096;
%    cdata4(idx,Element+1) = HIRSwords(frame,6)+4096;
%    cdata5(idx,Element+1) = HIRSwords(frame,7)+4096;
%    cdata6(idx,Element+1) = HIRSwords(frame,8)+4096;
%    cdata7(idx,Element+1) = HIRSwords(frame,9)+4096;
%    cdata8(idx,Element+1) = HIRSwords(frame,10)+4096;
%    cdata9(idx,Element+1) = HIRSwords(frame,11)+4096;
%    cdata10(idx,Element+1) = HIRSwords(frame,12)+4096;
%    cdata11(idx,Element+1) = HIRSwords(frame,13)+4096;
%    cdata12(idx,Element+1) = HIRSwords(frame,14)+4096;
%    cdata13(idx,Element+1) = HIRSwords(frame,15)+4096;
%    cdata14(idx,Element+1) = HIRSwords(frame,16)+4096;
%    cdata15(idx,Element+1) = HIRSwords(frame,17)+4096;
%    cdata16(idx,Element+1) = HIRSwords(frame,18)+4096;
%    cdata17(idx,Element+1) = HIRSwords(frame,19)+4096;
%    cdata18(idx,Element+1) = HIRSwords(frame,20)+4096;
%    cdata19(idx,Element+1) = HIRSwords(frame,21)+4096;
%    cdata20(idx,Element+1) = HIRSwords(frame,22)+4096;
   
   %if(mod(frame,62) == 0)
   %    idx = idx + 1;
   %end
   
   %if frameTime(frame)-frameTime(1)/12.6 == 0
   %   idx = idx + 1;
   %end
   
   %each scan takes 6.4 seconds
   idx = floor((frameTime(frame)-frameTime(1)+1.5)/(6.4))+1; 
   
   
   %idx = mod(frame,64)+1
   %idx = idx + 1;
   %for(idx = 1:20) %associate all 20 words with a specific element
   %   HIRSElements(frame,Element+1,idx) = HIRSwords(frame,idx+2);
   %end    
end


bright = 0;
contrast = 0;
%bright = 2.02^11;  %chan 22
%contrast = 2.19^6; %chan 22
%bright = 2.178^11;  %chan 21
%contrast = 2.235^6; %chan 21
%bright = 2.08^11;  %chan 20
%contrast = 2.233^6; %chan 20
%bright = 2.212^11;  %chan 19
%contrast = 2.211^6; %chan 19
%bright = 2.11^11;  %chan 18
%contrast = 2.211^6; %chan 18
%bright = 2.115^11;  %chan 17
%contrast = 2.175^6; %chan 17
%bright = 2.055^12;  %chan 16
%contrast = 2.17^6; %chan 16
%bright = 2^11;  %chan 15
%contrast = 2.25^5; %chan 15
%bright = 0;  %chan 14
%contrast = 2.05^6; %chan 14
%bright = 2^10;  %chan 13
%contrast = 2.25^5; %chan 13
%bright = 2.1^11;  %chan 12
%contrast = 2.11^6; %chan 12
%bright = 2^11;  %chan 11
%contrast = 2^5; %chan 11
%bright = 2^11; %chan 10
%contrast = 86;  %chan 10
%bright = 2.18^10; %chan 8
%contrast = 85;  %chan 8
%bright = 2.02^11; %channel 4
%contrast = 80; %channel 4

set(0,'DefaultFigureWindowStyle','docked')

for idx2 = 1:20
    %you should probably apply calibration data from the channels as well,
    %its in the user's guide
    
    HIRSstdDev = std(double(max(cdataTrusted(:,1:56,idx2)))); 
    HIRSmean = mean(double(max(cdataTrusted(:,1:56,idx2))));
    HIRSmax_nooutliers = max(abs((max(cdataTrusted(:,1:56,idx2)) > (2*HIRSstdDev+HIRSmean))-1) .* double(max(cdataTrusted(:,1:56,idx2))));
    
    figure(idx2);     
    %bright = max(max(cdata(:,1:53,1)));
    
    %contrast = max(max(cdataTrusted(:,1:56,idx2)))/1024
    %contrast = HIRSmax_nooutliers/255
    [row, col, v] = find(cdataTrusted(:,1:56,idx2));
    
    HIRSmin = mean(double(v))-2*std(double(v));
    HIRSmax = (HIRSmax_nooutliers-HIRSmin)/64;
    
    
    %%or if we know the calibration row...
    %HIRSmin =  mode (double(cdataTrusted(52,1:56,idx2)))
    %HIRSmin = HIRSmin-0.1*HIRSmin;
    %HIRSmax = mode (double(cdataTrusted(51,1:56,idx2))-HIRSmin)
    %HIRSmax = (HIRSmax+0.1*HIRSmax)/64;
    
    %*****pick one or the other depending on ascending or descending pass!
    %imgData = uint8(fliplr((cdata(:,1:maxElement,idx2)-HIRSmin)/(HIRSmax))); %DESCENDING
    imgData = uint8(flipud((cdata(:,1:maxElement,idx2)-HIRSmin)/(HIRSmax))); %ASCENDING
    
    image(imgData); 
    imwrite(imgData,colormap('bone'),['HIRS' num2str(idx2) '.png']);
    
    %image(flipud((cdata(:,1:56,idx2)-HIRSmin)/(HIRSmax))); 
    
    %image(flipud((cdataTrusted(:,1:56,idx2)-HIRSmin)/(HIRSmax)));
    colormap('bone');
    %colormap('summer');
    colorbar;
    pbaspect([size(cdata(:,:,1),2)*0.8 size(cdata(:,:,1),1) 1])
    saveas(idx2, ['HIRS_fig' num2str(idx2)], 'png')
end
%plot(cdata);

%% DCS - Pull out DCS-2 data in its own array
%18,19,
%24,25,
%28,29,
%32,33,
%40,41,
%44,45,
%52,53,
%56,57,
%60,61,
%64,65,
%68,69,
%72,73,
%76,77,
%86,87,
%90,91,
%94,95
clear DCSdata DCSTime;
DCSdata(:,1) = uint8(minorFrames(:,18+1)); DCSTime(:,1) = frameTime(:,18+1);
DCSdata(:,2) = (minorFrames(:,19+1)); DCSTime(:,2) = frameTime(:,19+1);

DCSdata(:,3) = (minorFrames(:,24+1)); DCSTime(:,3) = (frameTime(:,24+1));
DCSdata(:,4) = (minorFrames(:,25+1)); DCSTime(:,4) = (frameTime(:,25+1));

DCSdata(:,5) = (minorFrames(:,28+1)); DCSTime(:,5) = (frameTime(:,28+1));
DCSdata(:,6) = (minorFrames(:,29+1)); DCSTime(:,6) = (frameTime(:,29+1));

DCSdata(:,7) = (minorFrames(:,32+1)); DCSTime(:,7) = (frameTime(:,32+1));
DCSdata(:,8) = (minorFrames(:,33+1)); DCSTime(:,8) = (frameTime(:,33+1));

DCSdata(:,9) = (minorFrames(:,40+1)); DCSTime(:,9) = (frameTime(:,40+1));
DCSdata(:,10) = (minorFrames(:,41+1)); DCSTime(:,10) = (frameTime(:,41+1));

DCSdata(:,11) = (minorFrames(:,44+1)); DCSTime(:,11) = (frameTime(:,44+1));
DCSdata(:,12) = (minorFrames(:,45+1)); DCSTime(:,12) = (frameTime(:,45+1));

DCSdata(:,13) = (minorFrames(:,52+1)); DCSTime(:,13) = (frameTime(:,52+1));
DCSdata(:,14) = (minorFrames(:,53+1)); DCSTime(:,14) = (frameTime(:,53+1));

DCSdata(:,15) = (minorFrames(:,56+1)); DCSTime(:,15) = (frameTime(:,56+1));
DCSdata(:,16) = (minorFrames(:,57+1)); DCSTime(:,16) = (frameTime(:,57+1));

DCSdata(:,17) = (minorFrames(:,60+1)); DCSTime(:,17) = (frameTime(:,60+1));
DCSdata(:,18) = (minorFrames(:,61+1)); DCSTime(:,18) = (frameTime(:,61+1));

DCSdata(:,19) = (minorFrames(:,64+1)); DCSTime(:,19) = (frameTime(:,64+1));
DCSdata(:,20) = (minorFrames(:,65+1)); DCSTime(:,20) = (frameTime(:,65+1));

DCSdata(:,21) = (minorFrames(:,68+1)); DCSTime(:,21) = (frameTime(:,68+1));
DCSdata(:,22) = (minorFrames(:,69+1)); DCSTime(:,22) = (frameTime(:,69+1));

DCSdata(:,23) = (minorFrames(:,72+1)); DCSTime(:,23) = (frameTime(:,72+1));
DCSdata(:,24) = (minorFrames(:,73+1)); DCSTime(:,24) = (frameTime(:,73+1));

DCSdata(:,25) = (minorFrames(:,76+1)); DCSTime(:,25) = (frameTime(:,76+1));
DCSdata(:,26) = (minorFrames(:,77+1)); DCSTime(:,26) = (frameTime(:,77+1));

DCSdata(:,27) = (minorFrames(:,86+1)); DCSTime(:,27) = (frameTime(:,86+1));
DCSdata(:,28) = (minorFrames(:,87+1)); DCSTime(:,28) = (frameTime(:,87+1));

DCSdata(:,29) = (minorFrames(:,90+1)); DCSTime(:,29) = (frameTime(:,90+1));
DCSdata(:,30) = (minorFrames(:,91+1)); DCSTime(:,30) = (frameTime(:,91+1));

DCSdata(:,31) = (minorFrames(:,94+1)); DCSTime(:,31) = (frameTime(:,94+1));
DCSdata(:,32) = (minorFrames(:,95+1)); DCSTime(:,32) = (frameTime(:,95+1));

%% DCS - unroll DCS data into single stream and identify locations of embedded sync words
clear DCSstream packetstart DCSstreamTime;
%test(1)=uint16(0);
DCSstream(1)=uint8(0);
sFlag=0;
idx=1;
idx2=1;

if mode(spaceCraft) == 15
    byte2Max = 9;
else
    byte2Max = 8;
end

for frame=1:size(DCSdata,1)
    fprintf('\n');
    for byte=1:32        
          DCSstream(idx) = DCSdata(frame,byte);
          DCSstreamTime(idx) = DCSTime(frame,byte);
          fprintf('%0.2x ',DCSstream(idx));          
          if(sFlag > 0 && DCSdata(frame,byte) < byte2Max)
              packetstart(idx2) = sFlag;
              idx2 = idx2+1;
              sFlag = 0;
          else
              sFlag = 0;
          end
          
          if(DCSdata(frame,byte) == 214) 
              sFlag=idx-1;
          end
          
          idx = idx + 1;
    end
end
figure(6);
scatter(DCSstreamTime,DCSstream,2,'.');
axis tight;

%%%%%%%%%%%%%%%%%%%%%%%
% Table 1.2.2.6-1. DCS-2 System Characteristics.
% Parameter % Characteristic

% Minimum satellite elevation angle from platform % 5 degrees

% Number of platforms requiring location/velocity measurements visible in a 5 degree -visibility circle
% Capacity: 230

% Total number of such platforms over the globe
% Capacity: 4100

% Percentage of platforms with six good Doppler measurements per day
% 85%

% Platform transmission repetition period
% Approx. 60 sec

% Message length
% 0.3 to 0.9 sec

% Expected location accuracy
% 5 km to 8 km rms

% Expected velocity accuracy
% 1 to 1.6 mps

%%find packet header types 
%just an experiment, conclusion is there is only one sync word 'header'
%word
%clear headersx headersy;
%headersx(1)=0;
%headersy(1)=0;
%for idx=1:255
%    if(numel(strfind(DCSstream, [0 0 0 0 0 0 idx])) > 0)
%        %fprintf([num2str(idx) ' ' num2str(numel(strfind(DCSstream, [0 0 idx]))) '\n']);        
%        headersy(end+1)=numel(strfind(DCSstream, [0 0 idx]));
%        headersx(end+1)=idx;
%    end   
%end
%figure(7);
%stem(headersx, headersy);

%% DCS - Spit out raw hex values TO THE CONSOLE
clc;
for idx=1:numel(packetstart)-1     
    fprintf([num2str(idx) ' ' num2str(DCSstreamTime(packetstart(idx))) ' ' ]);          
    for byte=0:64             
        if(packetstart(idx)+byte == packetstart(idx+1))
            break;
        end                        
        fprintf(dec2hex(DCSstream(packetstart(idx)+byte+1),2));        
        fprintf(' ');
    end 
    fprintf('\n');     
end

%% DCS - Spit out raw hex DCS values TO A FILE

%byte 3 or 4
clc;
fid = fopen('DCS_RAW.txt', 'wt');
for idx=1:numel(packetstart)-1     
    fprintf(fid,'%04d %09.4f ',idx,DCSstreamTime(packetstart(idx)));
    fprintf('%04d %09.4f ',idx,DCSstreamTime(packetstart(idx)));
    for byte=0:64             
        if(packetstart(idx)+byte == packetstart(idx+1))
            break;
        end                        
        fprintf(dec2hex(DCSstream(packetstart(idx)+byte+1),2));        
        fprintf(fid, dec2hex(DCSstream(packetstart(idx)+byte+1),2));        
        fprintf(' ');
        fprintf(fid,' ');
    end 
    fprintf('\n');
    fprintf(fid,'\n');
end
fclose(fid);
% %%DCS -mess around with DCS ideas. Obsolete at this point. 
% %%%%%%%%%%%%%%%%%%%
% % 107 - value from counting signal, view surrounding bytes
% % 
% % 20406:  0    0    0  214    1   10   70  107   28  163    1    0    1    2    3    4
% % 20502:  0    0    0  214    7   91  198  107   44   54  222  192   76  246  186   64
% % 20582: 0    0    0  214    4  205  166  107   37  183  100   16  242   49  192   11
% %0 0 214 must be a sync word  
% 
% %sync+5 has these values. 
% % 6  38 70 102 134 166 198 230 
% % 00000110
% % 00100110
% % 01000110
% % 01100110
% % 10000110
% % 10100110
% % 11000110
% % 11100110
% 
% clear testy testx idx test;
% %testy(1,1)=uint32(0);
% %testx(1,1)=0;
% %idx2(8)=0;
% %idx3=1;
% %pattern = 'Ö(';
% %regexp(text, pattern, 'match')
% %packetstart = strfind(DCSstream, [214 ]); %Sync should be 0xD60 but this is hard in matlab because it crosses a byte boundary. 
% %packetstart = find(DCSstream == 214);
% %numel(packetstart);
% %NUT = 4;
% 
% %type 214
% for idx=1:numel(packetstart)-1 %s    
%     %packetstart(idx)+1 this is always 214
%     %packetstart(idx)+2 Upper nibble 0, lower nibble 0-7 (DCS channel?)
%     %packetstart(idx)+3 This appears to be a function of length, upper
%     %nibble:
%     %0=16 bytes
%     %3=20 bytes
%     %5=24 bytes
%     %6=28 bytes
%     %9=32 bytes
%     %A=36 bytes
%     %C=40 bytes
%     %F=44 bytes
%     
%     %packetstart(idx)+4 Lower 5 bits are the upper counter bits, others are possibly lower from previous byte?     
%     %packetstart(idx)+5 contains a counting up signal?
%     %packetstart(idx)+6 another count byte? (that would make this a 20-bit counter signal, average 104 counts/sec here, seems OK as time axis is approximate)
%     %packetstart(idx)+7 contains embedded trends (ID byte 1)
%     %packetstart(idx)+8 contains random data (ID byte 2)
%     %packetstart(idx)+9 contains 0-16 in the upper nibble with 0's in the
%     %lower nibble (ID byte 3)
%     %packetstart(idx)+10 contains almost evenly spaced values 
%     %packetstart(idx)+11 contains multiple data sets with most between 0-64    
%     %packetstart(idx)+12 contains random data 
%     %packetstart(idx)+13 contains random data
%     %packetstart(idx)+14 contains random data
%     %packetstart(idx)+15 contains random data        
%     
%     %%%fprintf([num2str(idx) ' ' num2str(DCSstreamTime(packetstart(idx))) ' ' ]);    
%     %DCSPacket(byte+1,idx)=DCSstream(packetstart(idx)+byte+1);
%     
%     %bitshift(128,8,16)
%     %bitshift(bitand(127,bin2dec('11111')),8,16)
%     %bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))
%     %bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),8,32)
%     
%     DCSPacketTime = (bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + bitshift(double(DCSstream(packetstart(idx)+5)),8,32) + double(DCSstream(packetstart(idx)+6)));
%     fprintf(num2str(DCSPacketTime));
%     fprintf(' ');
%     fprintf(dec2bin(DCSPacketTime));
%     fprintf('\n');
%     DCSPacketTime = DCSPacketTime/100.0;
%     
% %     for byte=0:64             
% %         if(packetstart(idx)+byte == packetstart(idx+1))
% %             break;
% %         end
% %         
% %         if(byte == 3-1)
% %            testy2(idx)= DCSstream(packetstart(idx)+byte+1)
% %         elseif(byte == 4-1)
% %            testx(idx)=DCSstreamTime(packetstart(idx)+3);
% %            testy2(idx)= bitand(DCSstream(packetstart(idx)+byte+1),bin2dec('11111'));
% %            testy(idx)= bitand(DCSstream(packetstart(idx)+byte+1),bin2dec('11111'));
% %            %testy(idx)=0;
% %         elseif byte == 5-1
% %            testy(idx) = bitshift(testy(idx),8)+uint32(DCSstream(packetstart(idx)+byte+1));
% %         elseif byte == 6-1
% %             testy(idx) = bitshift(testy(idx),8)+uint32(DCSstream(packetstart(idx)+byte+1));
% %             %testy(idx) = bitand(uint16(DCSstream(packetstart(idx)+byte+1)),bin2dec('1111'));
% %         end
% %         
% %         %%%fprintf(dec2hex(DCSstream(packetstart(idx)+byte+1),2));
% %         %DCSPacket(byte+1,idx)=DCSstream(packetstart(idx)+byte+1);
% %         %fprintf(num2str(DCSstream(packetstart(idx)+byte+1)));
% %         %%%fprintf(' ');
% %     end 
%     %%%fprintf('\n');     
% end
% %plot(test(:,:),'.-');
% %chan=2;
% %stairs(testx,testy,'.-');
% %stairs(test,'.-');
% %figure(8);
% %plot(testx,testy,'-x');
% %plot(test(:,1),test(:,2),'.-');
% %axis([1 max(test(:,1)) 0 255]);

%% DCS - Write to files based on TXID
%%WARNING -- THIS WILL WRITE A BUTTLOAD OF TINY FILES!
%%PLEASE SET YOUR CURRENT FOLDER ACCORDINGLY

%Beware! ADCS spacecraft (IE NOAA19) have an additional byte and encodes
%lengthes differently than DCS2!

%ADCS - byte 2 is length
%1-16
%2-20
%3-24
%4-28
%5-32
%6-36
%7-40
dontWriteToFiles = 0;
dontWriteSummaryFile = 0;

TXIDList = containers.Map; %Hash list that contains all transmitters


for idx=1:numel(packetstart)-1     
    TXID = [dec2hex(DCSstream(packetstart(idx)-1+8),2) dec2hex(DCSstream(packetstart(idx)-1+9),2) dec2hex(DCSstream(packetstart(idx)-1+10),2) dec2hex(DCSstream(packetstart(idx)-1+11),2)];
    
    if(packetstart(idx)-1+11 >= packetstart(idx+1))
       continue;
    end
    
    %DCSPacketTime = (bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + double(bitshift(DCSstream(packetstart(idx)+5),8,32)) + double(DCSstream(packetstart(idx)+6))) / 100.0;
    
    fprintf([TXID ',']);
    fprintf([dec2hex(bitshift(DCSstream(packetstart(idx)+3),-4)) '\n']);
    
    switch(dec2hex(bitshift(DCSstream(packetstart(idx)+3),-4)))
        case '0'
            numbytes=16;
        case '3'
            numbytes=20;
        case '5'
            numbytes=24;
        case '6'
            numbytes=28;
        case '9'
            numbytes=32;
        case 'A'
            numbytes=36;
        case 'C'
            numbytes=40;
        case 'F'
            numbytes=44;
        otherwise 
            numbytes=44;
    end
    
    if dontWriteToFiles == 0
        fid = fopen([TXID '.txt'],'at');
        fprintf(fid,'%04d %09.4f ',idx,DCSstreamTime(packetstart(idx)));
    end
    
    %print ID and session time
    %fprintf(fid,[num2str(idx) ' ' num2str(DCSstreamTime(packetstart(idx))) ' ' ]);
    %fprintf([num2str(idx) ' ' num2str(DCSstreamTime(packetstart(idx))) ' ' ]);          
    fprintf('%04d %09.4f ',idx,DCSstreamTime(packetstart(idx)));    
    
    for byte=0:numbytes-1
        if(packetstart(idx)+byte == packetstart(idx+1))
            
            break;
        end                        
        fprintf(dec2hex(DCSstream(packetstart(idx)+byte+1),2));
        fprintf(' ');
        
        if dontWriteToFiles == 0
            fprintf(fid, dec2hex(DCSstream(packetstart(idx)+byte+1),2));        
            fprintf(fid,' ');
        end
        
    end
    %fprintf(num2str((bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + bitshift(double(DCSstream(packetstart(idx)+5)),8,32) + double(DCSstream(packetstart(idx)+6))) / 128.0));
    %fprintf(num2str((bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + bitshift(double(DCSstream(packetstart(idx)+5)),8,32) + double(DCSstream(packetstart(idx)+6))) / 100.0));
    
    %fprintf(num2str(0.9549*(bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,'uint32') + bitshift(double(DCSstream(packetstart(idx)+5)),8,'uint32') + double(DCSstream(packetstart(idx)+6))) / 100.0));
    %fprintf(num2str(0.9549*(bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + bitshift(double(DCSstream(packetstart(idx)+5)),8,32) + double(DCSstream(packetstart(idx)+6))) / 100.0));
    %fprintf(' ');
    
    %fprintf(fid,num2str((bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + bitshift(double(DCSstream(packetstart(idx)+5)),8,32) + double(DCSstream(packetstart(idx)+6))) / 128.0));
    %fprintf(fid,num2str((bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + bitshift(double(DCSstream(packetstart(idx)+5)),8,32) + double(DCSstream(packetstart(idx)+6))) / 100.0));
    %fprintf(fid,num2str(0.9549*(bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,'uint32') + bitshift(double(DCSstream(packetstart(idx)+5)),8,'uint32') + double(DCSstream(packetstart(idx)+6))) / 100.0));
    
    %Different versions of matlab handle the third parameter differently, older want integer, newer want uint32 string
    if dontWriteToFiles == 0
       fprintf(fid,'%011.5f ',0.9549*(bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,'uint32') + bitshift(double(DCSstream(packetstart(idx)+5)),8,'uint32') + double(DCSstream(packetstart(idx)+6))) / 100.0);
       %fprintf(fid,'%011.5f ',0.9549*(bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + bitshift(double(DCSstream(packetstart(idx)+5)),8,32) + double(DCSstream(packetstart(idx)+6))) / 100.0);
    end
    
    fprintf('%011.5f ',0.9549*(bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,'uint32') + bitshift(double(DCSstream(packetstart(idx)+5)),8,'uint32') + double(DCSstream(packetstart(idx)+6))) / 100.0);
    %fprintf('%011.5f ',0.9549*(bitshift(double(bitand(DCSstream(packetstart(idx)+4),bin2dec('11111'))),16,32) + bitshift(double(DCSstream(packetstart(idx)+5)),8,32) + double(DCSstream(packetstart(idx)+6))) / 100.0);
        
    %parity bit should ALWAYS make this even, if odd, then problem. Mark it as such. 
    if mod(numel(find(dec2bin(hex2dec([dec2hex(DCSstream(packetstart(idx)+numbytes-2),2) dec2hex(DCSstream(packetstart(idx)+numbytes-1),2) dec2hex(DCSstream(packetstart(idx)+numbytes),2)]))=='1')),2) ~= 0
        fprintf('*');
        if dontWriteToFiles == 0
            fprintf(fid,'*');
        end
    end
    
    %Use all 23 bits including 5 fractional bits (last one is parity, dealt with it above)
    %fprintf(num2str((bitshift(hex2dec([dec2hex(DCSstream(packetstart(idx)+numbytes-2),2) dec2hex(DCSstream(packetstart(idx)+numbytes-1),2) dec2hex(DCSstream(packetstart(idx)+numbytes),2)]),-1)-2^22)/32.0));
    %fprintf(fid,num2str((bitshift(hex2dec([dec2hex(DCSstream(packetstart(idx)+numbytes-2),2) dec2hex(DCSstream(packetstart(idx)+numbytes-1),2) dec2hex(DCSstream(packetstart(idx)+numbytes),2)]),-1)-2^22)/32.0));
    fprintf('%012.5f ',(bitshift(hex2dec([dec2hex(DCSstream(packetstart(idx)+numbytes-2),2) dec2hex(DCSstream(packetstart(idx)+numbytes-1),2) dec2hex(DCSstream(packetstart(idx)+numbytes),2)]),-1)-2^22)/32.0);
    fprintf('\n');
    
    if dontWriteToFiles == 0
        fprintf(fid,'%012.5f ',(bitshift(hex2dec([dec2hex(DCSstream(packetstart(idx)+numbytes-2),2) dec2hex(DCSstream(packetstart(idx)+numbytes-1),2) dec2hex(DCSstream(packetstart(idx)+numbytes),2)]),-1)-2^22)/32.0);
        fprintf(fid,'\n');
        fclose(fid);
    end
    
    %doesn't grab last 6 bits
    %fprintf(num2str(bitshift(hex2dec([dec2hex(DCSstream(packetstart(idx)+numbytes-2),2) dec2hex(DCSstream(packetstart(idx)+numbytes-1),2) dec2hex(DCSstream(packetstart(idx)+numbytes),2)]),-6)-2^17));
    %fprintf(fid,num2str(bitshift(hex2dec([dec2hex(DCSstream(packetstart(idx)+numbytes-2),2) dec2hex(DCSstream(packetstart(idx)+numbytes-1),2) dec2hex(DCSstream(packetstart(idx)+numbytes),2)]),-6)-2^17));
        
    if TXIDList.isKey(TXID)
        TXIDList(TXID) = TXIDList(TXID) + 1;
    else
        TXIDList(TXID) = 1;
    end  
end

TXIDKeys = TXIDList.keys;
TXIDcount = cell2mat(TXIDList.values);

[TXIDcount, sortIndex] = sort(TXIDcount,'descend');
TXIDKeys = TXIDKeys(sortIndex);


if dontWriteSummaryFile == 0
   fid = fopen('DCS_summary.txt', 'wt');
end

for loop=1:numel(TXIDKeys)
    fprintf('%s %d \n',cell2mat(TXIDKeys(loop)), TXIDcount(loop));
    if dontWriteSummaryFile == 0      
       fprintf(fid, '%s %d \n',cell2mat(TXIDKeys(loop)), TXIDcount(loop));
    end
end

if dontWriteSummaryFile == 0
    fclose(fid);
end


%% SEM - Pull SEM bytes from format into two streams
%(one stream for each embedded byte)
minorFrameID = bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6));
SEMdata(:,1) = 255 - (minorFrames(:,20+1)); %SEM data appears to be inverted! *is it possible that MSB and LSB are interchanged?
SEMdata(:,2) = 255 - (minorFrames(:,21+1)); %SEM data appears to be inverted! *Double check this (parity calcs showed this to be an issue)

%% SEM - Pull out SEM data MEPED "Medium Energy Proton and Electron Detector"
%MEPED Digital A data consists of six directional proton measurements and three directional electron measurements for each of two directions of incidence (0 and 90 degrees) and four omni-directional proton measurements. All but the two highest energy omni-directional proton measurements are read out every two seconds. The two highest energy omnidirectional proton measurements are read out every four seconds. The MEPED Digital A data and readout rates are summarized in Table 4.3.4.2-2.
clear MEPED_0P1 MEPED_0P2 MEPED_0P3 MEPED_0P4 MEPED_0P5 MEPED_0P6;
clear MEPED_0E1 MEPED_0E2 MEPED_0E3 MEPED_9E1 MEPED_9E2 MEPED_9E3
clear MEPED_9P1 MEPED_9P2 MEPED_9P3 MEPED_9P4 MEPED_9P5 MEPED_9P6;
clear MEPED_P6 MEPED_P7 MEPED_P8 MEPED_P9

MEPED_0P1(1)=0; MEPED_0P2(1)=0; MEPED_0P3(1)=0; MEPED_0P4(1)=0; MEPED_0P5(1)=0; MEPED_0P6(1)=0;
MEPED_9P1(1)=0; MEPED_9P2(1)=0; MEPED_9P3(1)=0; MEPED_9P4(1)=0; MEPED_9P5(1)=0; MEPED_9P6(1)=0;
MEPED_0E1(1)=0; MEPED_0E2(1)=0; MEPED_0E3(1)=0; MEPED_9E1(1)=0; MEPED_9E2(1)=0; MEPED_9E3(1)=0;
MEPED_P6(1)=0; MEPED_P7(1)=0; MEPED_P8(1)=0; MEPED_P9(1)=0;


for frame=1:numel(minorFrameID)
    %if(minorFrameID(frame) == 0)
    %    day=bitshift(minorFrames(frame,8+1),1)+bitshift(bitor(minorFrames(frame,9+1),128),-7); %should be 241 for loproto2 recording        
        %fprintf([num2str(day) ' ' num2str(dec2hex(minorFrames(frame,8+1))) ' ' num2str(dec2hex(minorFrames(frame,9+1))) '\n\n']);
    %    fprintf([num2str(day) '\n\n']);
        
    %end
    if(mod(minorFrameID(frame),20) == 0)
        MEPED_0P1(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-1,20) == 0)
        MEPED_0P2(end+1) = SEMdata(frame,1);
        MEPED_0P3(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-2,20) == 0)
        MEPED_0P4(end+1) = SEMdata(frame,1);
        MEPED_0P5(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-3,20) == 0)
        MEPED_0P6(end+1) = SEMdata(frame,1);
        MEPED_0E1(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-4,20) == 0)
        MEPED_0E2(end+1) = SEMdata(frame,1);
        MEPED_0E3(end+1) = SEMdata(frame,2);    
    elseif(mod(minorFrameID(frame)-5,20) == 0)
        MEPED_9P1(end+1) = SEMdata(frame,1);
        MEPED_9P2(end+1) = SEMdata(frame,2);    
    elseif(mod(minorFrameID(frame)-6,20) == 0)
        MEPED_9P3(end+1) = SEMdata(frame,1);
        MEPED_9P4(end+1) = SEMdata(frame,2);        
    elseif(mod(minorFrameID(frame)-7,20) == 0)
        MEPED_9P5(end+1) = SEMdata(frame,1);
        MEPED_9P6(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-8,20) == 0)
        MEPED_9E1(end+1) = SEMdata(frame,1);
        MEPED_9E2(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-9,20) == 0)
        MEPED_9E3(end+1) = SEMdata(frame,1);
        MEPED_P6(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-10,20) == 0)
        MEPED_P7(end+1) = SEMdata(frame,1);        
    end
    if(mod(minorFrameID(frame)-10,40) == 0)        
        MEPED_P8(end+1) = SEMdata(frame,2);
    end
    if(mod(minorFrameID(frame)-30,40) == 0)        
        MEPED_P9(end+1) = SEMdata(frame,2);
    end
end
%% SEM - Filter MEPED data
%remove bits that change too much over a single period (fail a lookhead and
%behind test)
FilterThreshold = 20;

for idx=2:numel(MEPED_0E1)-1
   if(abs(MEPED_0E1(idx-1) - MEPED_0E1(idx)) > FilterThreshold && abs(MEPED_0E1(idx+1) - MEPED_0E1(idx)) > FilterThreshold)
      MEPED_0E1(idx) = 0;
   end
end

for idx=2:numel(MEPED_0E2)-1
   if(abs(MEPED_0E2(idx-1) - MEPED_0E2(idx)) > FilterThreshold && abs(MEPED_0E2(idx+1) - MEPED_0E2(idx)) > FilterThreshold)
      MEPED_0E2(idx) = 0;
   end
end

for idx=2:numel(MEPED_0E3)-1
  if(abs(MEPED_0E3(idx-1) - MEPED_0E3(idx)) > FilterThreshold && abs(MEPED_0E3(idx+1) - MEPED_0E3(idx)) > FilterThreshold)
      MEPED_0E3(idx) = 0;
   end
end

for idx=2:numel(MEPED_9E1)-1
   if(abs(MEPED_9E1(idx-1) - MEPED_9E1(idx)) > FilterThreshold && abs(MEPED_9E1(idx+1) - MEPED_9E1(idx)) > FilterThreshold)
      MEPED_9E1(idx) = 0;
   end
end

for idx=2:numel(MEPED_9E2)-1
   if(abs(MEPED_9E2(idx-1) - MEPED_9E2(idx)) > FilterThreshold && abs(MEPED_9E2(idx+1) - MEPED_9E2(idx)) > FilterThreshold)
      MEPED_9E2(idx) = 0;
   end
end

for idx=2:numel(MEPED_9E3)-1
   if(abs(MEPED_9E3(idx-1) - MEPED_9E3(idx)) > FilterThreshold && abs(MEPED_9E3(idx+1) - MEPED_9E3(idx)) > FilterThreshold)
      MEPED_9E3(idx) = 0;
   end
end

for idx=2:numel(MEPED_0P1)-1
   if(abs(MEPED_0P1(idx-1) - MEPED_0P1(idx)) > FilterThreshold && abs(MEPED_0P1(idx+1) - MEPED_0P1(idx)) > FilterThreshold)
      MEPED_0P1(idx) = 0;
   end
end

for idx=2:numel(MEPED_0P2)-1
   if(abs(MEPED_0P2(idx-1) - MEPED_0P2(idx)) > FilterThreshold && abs(MEPED_0P2(idx+1) - MEPED_0P2(idx)) > FilterThreshold)
      MEPED_0P2(idx) = 0;
   end
end

for idx=2:numel(MEPED_0P3)-1
  if(abs(MEPED_0P3(idx-1) - MEPED_0P3(idx)) > FilterThreshold && abs(MEPED_0P3(idx+1) - MEPED_0P3(idx)) > FilterThreshold)
      MEPED_0P3(idx) = 0;
   end
end

for idx=2:numel(MEPED_0P4)-1
   if(abs(MEPED_0P4(idx-1) - MEPED_0P4(idx)) > FilterThreshold && abs(MEPED_0P4(idx+1) - MEPED_0P4(idx)) > FilterThreshold)
      MEPED_0P4(idx) = 0;
   end
end

for idx=2:numel(MEPED_0P5)-1
   if(abs(MEPED_0P5(idx-1) - MEPED_0P5(idx)) > FilterThreshold && abs(MEPED_0P5(idx+1) - MEPED_0P5(idx)) > FilterThreshold)
      MEPED_0P5(idx) = 0;
   end
end

for idx=2:numel(MEPED_0P6)-1
   if(abs(MEPED_0P6(idx-1) - MEPED_0P6(idx)) > FilterThreshold && abs(MEPED_0P6(idx+1) - MEPED_0P6(idx)) > FilterThreshold)
      MEPED_0P6(idx) = 0;
   end
end

for idx=2:numel(MEPED_9P1)-1
   if(abs(MEPED_9P1(idx-1) - MEPED_9P1(idx)) > FilterThreshold && abs(MEPED_9P1(idx+1) - MEPED_9P1(idx)) > FilterThreshold)
      MEPED_9P1(idx) = 0;
   end
end

for idx=2:numel(MEPED_9P2)-1
   if(abs(MEPED_9P2(idx-1) - MEPED_9P2(idx)) > FilterThreshold && abs(MEPED_9P2(idx+1) - MEPED_9P2(idx)) > FilterThreshold)
      MEPED_9P2(idx) = 0;
   end
end

for idx=2:numel(MEPED_9P3)-1
  if(abs(MEPED_9P3(idx-1) - MEPED_9P3(idx)) > FilterThreshold && abs(MEPED_9P3(idx+1) - MEPED_9P3(idx)) > FilterThreshold)
      MEPED_9P3(idx) = 0;
   end
end

for idx=2:numel(MEPED_9P4)-1
   if(abs(MEPED_9P4(idx-1) - MEPED_9P4(idx)) > FilterThreshold && abs(MEPED_9P4(idx+1) - MEPED_9P4(idx)) > FilterThreshold)
      MEPED_9P4(idx) = 0;
   end
end

for idx=2:numel(MEPED_9P5)-1
   if(abs(MEPED_9P5(idx-1) - MEPED_9P5(idx)) > FilterThreshold && abs(MEPED_9P5(idx+1) - MEPED_9P5(idx)) > FilterThreshold)
      MEPED_9P5(idx) = 0;
   end
end

for idx=2:numel(MEPED_9P6)-1
   if(abs(MEPED_9P6(idx-1) - MEPED_9P6(idx)) > FilterThreshold && abs(MEPED_9P6(idx+1) - MEPED_9P6(idx)) > FilterThreshold)
      MEPED_9P6(idx) = 0;
   end
end

%% SEM - Plot MEPED data
figure(9);

subplot(3,2,1);
plot(1:numel(MEPED_0E1),MEPED_0E1,1:numel(MEPED_0E2),MEPED_0E2,1:numel(MEPED_0E3),MEPED_0E3);
title('0 Degree Electron Count');
h_legend=legend('>=30 keV Electrons/2sec','>=100 keV Electrons/2sec','>=300 keV Electrons/2sec');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(MEPED_P6) -0.5 256]);

subplot(3,2,2);
plot(1:numel(MEPED_9E1),MEPED_9E1,1:numel(MEPED_9E2),MEPED_9E2,1:numel(MEPED_9E3),MEPED_9E3);
title('90 Degree Electron Count');
h_legend=legend('>=30 keV Electrons/2sec','>=100 keV Electrons/2sec','>=300 keV Electrons/2sec');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(MEPED_P6) -0.5 256]);

subplot(3,2,3);
plot(1:numel(MEPED_0P1),MEPED_0P1,1:numel(MEPED_0P2),MEPED_0P2,1:numel(MEPED_0P3),MEPED_0P3,1:numel(MEPED_0P4),MEPED_0P4,1:numel(MEPED_0P5),MEPED_0P5,1:numel(MEPED_0P6),MEPED_0P6);
title('0 Degree Proton Count');
h_legend=legend('30-80keV Protons/2sec','80-250keV Protons/2sec','250-800keV Protons/2sec','800-2.5k keV Protons/2sec','2.5k-7k keV Protons/2sec','>7k keV Protons/2sec');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(MEPED_P6) -0.5 256]);

subplot(3,2,4);
plot(1:numel(MEPED_9P1),MEPED_9P1,1:numel(MEPED_9P2),MEPED_9P2,1:numel(MEPED_9P3),MEPED_9P3,1:numel(MEPED_9P4),MEPED_9P4,1:numel(MEPED_9P5),MEPED_9P5,1:numel(MEPED_9P6),MEPED_9P6);
title('90 Degree Proton Count');
h_legend=legend('30-80keV Protons/2sec','80-250keV Protons/2sec','250-800keV Protons/2sec','800-2.5k keV Protons/2sec','2.5k-7k keV Protons/2sec','>7k keV Protons/2sec');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(MEPED_P6) -0.5 256]);

subplot(3,1,3);
plot(1:numel(MEPED_P6),MEPED_P6,1:numel(MEPED_P7),MEPED_P7,1:2:numel(MEPED_P8)*2,MEPED_P8,1:2:numel(MEPED_P9)*2,MEPED_P9);
title('Omni Proton Count');
h_legend=legend('>=16 MeV Protons/2sec','>=35 Mev Protons/2sec','>=70 Mev Protons/2sec','>=140 Mev Protons/2sec');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(MEPED_P6) -0.5 256]);

suptitle('SEM: Multiple Energy Proton/Electron Detector');

%% SEM - Pull out SEM data TED "Total Energy Detector"
%TED Digital A data consists of a 0.05 to 1 keV partial energy flux measurement, a 1 to 20 keV partial energy flux measurement, maximum differential energy fluxes, four-point differential energy spectra and background measurements for electrons and protons, each at two angles of incidence (0 and 30 degrees).

clear TED_0EFL TED_0PFL TED_3EFL TED_3PFL;
clear TED_0EFH TED_0PFH TED_3EFH TED_3PFH;
clear TED_0DEM TED_0DPM TED_3DEM TED_3DPM;
clear TED_0EM TED_0PM TED_3EM TED_3PM;
clear TED_0DE1 TED_0DE2 TED_0DE3 TED_0DE4 TED_3DE1 TED_3DE2 TED_3DE3 TED_3DE4 TED_0DP1 TED_0DP2 TED_0DP3 TED_0DP4 TED_3DP1 TED_3DP2 TED_3DP3 TED_3DP4;
clear TED_0EBKH TED_0EBKL TED_0PBKH TED_0PBKL TED_3PBKH TED_3PBKL;

TED_0EFL(1)=0;TED_0PFL(1)=0;TED_3EFL(1)=0;TED_3PFL(1)=0;
TED_0EFH(1)=0;TED_0PFH(1)=0;TED_3EFH(1)=0;TED_3PFH(1)=0;
TED_0DEM(1)=0;TED_0DPM(1)=0;TED_3DEM(1)=0;TED_3DPM(1)=0;
TED_0EM(1)=0;TED_0PM(1)=0;TED_3EM(1)=0;TED_3PM(1)=0;
TED_0DE1(1)=0; TED_0DE2(1)=0; TED_0DE3(1)=0; TED_0DE4(1)=0; TED_3DE1(1)=0; TED_3DE2(1)=0; TED_3DE3(1)=0; TED_3DE4(1)=0; TED_0DP1(1)=0; TED_0DP2(1)=0; TED_0DP3(1)=0; TED_0DP4(1)=0; TED_3DP1(1)=0; TED_3DP2(1)=0; TED_3DP3(1)=0; TED_3DP4(1)=0;
TED_0EBKH(1)=0;TED_0EBKL(1)=0;TED_0PBKH(1)=0;TED_0PBKL(1)=0;TED_3PBKH(1)=0;TED_3PBKL(1)=0;

%minorFrameID = bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6));
%DCSdata(:,1) = (minorFrames(:,20+1));
%DCSdata(:,2) = (minorFrames(:,21+1));
for frame=1:numel(minorFrameID)
    %if(minorFrameID(frame) == 0)
    %    day=bitshift(minorFrames(frame,8+1),1)+bitor(minorFrames(frame,9+1),128)
    %end
    if(mod(minorFrameID(frame)-13,20) == 0)
        TED_0EFL(end+1) = SEMdata(frame,1);
        TED_3EFL(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-14,20) == 0)
        TED_0PFL(end+1) = SEMdata(frame,1);
        TED_3PFL(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-15,20) == 0)
        TED_0EFH(end+1) = SEMdata(frame,1);
        TED_3EFH(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-16,20) == 0)
        TED_0PFH(end+1) = SEMdata(frame,1);
        TED_3PFH(end+1) = SEMdata(frame,2);        
    elseif(mod(minorFrameID(frame)-17,20) == 0)
        TED_0EM(end+1) = bitshift(bitand(SEMdata(frame,1),128+64+32+16),-4);
        TED_0PM(end+1) = bitand(SEMdata(frame,1),1+2+4+8);
        TED_0DEM(end+1) = SEMdata(frame,2);        
    elseif(mod(minorFrameID(frame)-18,20) == 0)
        TED_0DPM(end+1) = SEMdata(frame,1);
        TED_3EM(end+1) = bitshift(bitand(SEMdata(frame,2),128+64+32+16),-4);
        TED_3PM(end+1) = bitand(SEMdata(frame,2),1+2+4+8);
    elseif(mod(minorFrameID(frame)-19,20) == 0)
        TED_3DEM(end+1) = SEMdata(frame,1);
        TED_3DPM(end+1) = SEMdata(frame,2);        
    elseif(mod(minorFrameID(frame)-11,80) == 0)
        TED_0DE1(end+1) = SEMdata(frame,1);
        TED_0DE2(end+1) = SEMdata(frame,2);        
    elseif(mod(minorFrameID(frame)-31,80) == 0)
        TED_3DE1(end+1) = SEMdata(frame,1);
        TED_3DE2(end+1) = SEMdata(frame,2);        
    elseif(minorFrameID(frame) == 51 || minorFrameID(frame) == 131 || minorFrameID(frame) == 211)
        TED_0DP1(end+1) = SEMdata(frame,1);
        TED_0DP2(end+1) = SEMdata(frame,2);        
    elseif(minorFrameID(frame) == 71 || minorFrameID(frame) == 151 || minorFrameID(frame) == 231)
        TED_3DP1(end+1) = SEMdata(frame,1);
        TED_3DP2(end+1) = SEMdata(frame,2);        
    elseif(minorFrameID(frame) == 291)
        TED_0EBKL(end+1) = SEMdata(frame,1);
        TED_0EBKH(end+1) = SEMdata(frame,2);
    elseif(minorFrameID(frame) == 311)        
        %dex2hex(SEMdata(frame,1)) %should be 0xF3, data seems inverted?        
        %frame
        TED_3PBKL(end+1) = SEMdata(frame,2);
    elseif(mod(minorFrameID(frame)-12,80) == 0)
        TED_0DE3(end+1) = SEMdata(frame,1);
        TED_0DE4(end+1) = SEMdata(frame,2);        
    elseif(mod(minorFrameID(frame)-32,80) == 0)
        TED_3DE3(end+1) = SEMdata(frame,1);
        TED_3DE4(end+1) = SEMdata(frame,2);
    elseif(minorFrameID(frame) == 52 || minorFrameID(frame) == 132 || minorFrameID(frame) == 212)
        TED_0DP3(end+1) = SEMdata(frame,1);
        TED_0DP4(end+1) = SEMdata(frame,2);
    elseif(minorFrameID(frame) == 72 || minorFrameID(frame) == 152 || minorFrameID(frame) == 232)
        TED_3DP3(end+1) = SEMdata(frame,1);
        TED_3DP4(end+1) = SEMdata(frame,2);        
   elseif(minorFrameID(frame) == 292)
        TED_0PBKL(end+1) = SEMdata(frame,1);
        TED_0PBKH(end+1) = SEMdata(frame,2);    
    elseif(minorFrameID(frame) == 312)        
        %dex2hex(SEMdata(frame,1)) %should be 0x50, data is clearly inverted
        TED_3PBKH(end+1) = SEMdata(frame,2);
    end
    
end

%% SEM - Filter TED data
%remove bits that change too much over a single period (fail a lookhead and
%behind test)
FilterThreshold = 20;

for idx=2:numel(TED_0EFL)-1
   if(abs(TED_0EFL(idx-1) - TED_0EFL(idx)) > FilterThreshold && abs(TED_0EFL(idx+1) - TED_0EFL(idx)) > FilterThreshold)
      TED_0EFL(idx) = 0;
   end
end

for idx=2:numel(TED_0PFL)-1
   if(abs(TED_0PFL(idx-1) - TED_0PFL(idx)) > FilterThreshold && abs(TED_0PFL(idx+1) - TED_0PFL(idx)) > FilterThreshold)
      TED_0PFL(idx) = 0;
   end
end

for idx=2:numel(TED_3EFL)-1
  if(abs(TED_3EFL(idx-1) - TED_3EFL(idx)) > FilterThreshold && abs(TED_3EFL(idx+1) - TED_3EFL(idx)) > FilterThreshold)
      TED_3EFL(idx) = 0;
   end
end

for idx=2:numel(TED_3PFL)-1
   if(abs(TED_3PFL(idx-1) - TED_3PFL(idx)) > FilterThreshold && abs(TED_3PFL(idx+1) - TED_3PFL(idx)) > FilterThreshold)
      TED_3PFL(idx) = 0;
   end
end

for idx=2:numel(TED_0EFH)-1
   if(abs(TED_0EFH(idx-1) - TED_0EFH(idx)) > FilterThreshold && abs(TED_0EFH(idx+1) - TED_0EFH(idx)) > FilterThreshold)
      TED_0EFH(idx) = 0;
   end
end

for idx=2:numel(TED_0PFH)-1
   if(abs(TED_0PFH(idx-1) - TED_0PFH(idx)) > FilterThreshold && abs(TED_0PFH(idx+1) - TED_0PFH(idx)) > FilterThreshold)
      TED_0PFH(idx) = 0;
   end
end

for idx=2:numel(TED_3EFH)-1
  if(abs(TED_3EFH(idx-1) - TED_3EFH(idx)) > FilterThreshold && abs(TED_3EFH(idx+1) - TED_3EFH(idx)) > FilterThreshold)
      TED_3EFH(idx) = 0;
   end
end

for idx=2:numel(TED_3PFH)-1
   if(abs(TED_3PFH(idx-1) - TED_3PFH(idx)) > FilterThreshold && abs(TED_3PFH(idx+1) - TED_3PFH(idx)) > FilterThreshold)
      TED_3PFH(idx) = 0;
   end
end

for idx=2:numel(TED_0DEM)-1
   if(abs(TED_0DEM(idx-1) - TED_0DEM(idx)) > FilterThreshold && abs(TED_0DEM(idx+1) - TED_0DEM(idx)) > FilterThreshold)
      TED_0DEM(idx) = 0;
   end
end

for idx=2:numel(TED_0DPM)-1
   if(abs(TED_0DPM(idx-1) - TED_0DPM(idx)) > FilterThreshold && abs(TED_0DPM(idx+1) - TED_0DPM(idx)) > FilterThreshold)
      TED_0DPM(idx) = 0;
   end
end

for idx=2:numel(TED_3DEM)-1
  if(abs(TED_3DEM(idx-1) - TED_3DEM(idx)) > FilterThreshold && abs(TED_3DEM(idx+1) - TED_3DEM(idx)) > FilterThreshold)
      TED_3DEM(idx) = 0;
   end
end

for idx=2:numel(TED_3DPM)-1
   if(abs(TED_3DPM(idx-1) - TED_3DPM(idx)) > FilterThreshold && abs(TED_3DPM(idx+1) - TED_3DPM(idx)) > FilterThreshold)
      TED_3DPM(idx) = 0;
   end
end

%% SEM - Plot TED data
%0.05-1 keV Partial Energy Flux
%TED_0EFL TED_3EFL 
%TED_0PFL TED_3PFL

%2-10 keV Partial Energy Flux
%TED_0EFH TED_3EFH 
%TED_0PFH TED_3PFH

%Maximum Differential Energy Flux
%TED_0DEM  TED_3DEM 
%TED_0DPM TED_3DPM

%Energy of Maximum Differential Energy Flux
%TED_0EM TED_3EM 
%TED_0PM TED_3PM

%Four Point Energy/Flux Spectrum
%TED_0DE1 TED_0DE2 TED_0DE3 TED_0DE4 TED_3DE1 TED_3DE2 TED_3DE3 TED_3DE4 
%TED_0DP1 TED_0DP2 TED_0DP3 TED_0DP4 TED_3DP1 TED_3DP2 TED_3DP3 TED_3DP4

%Background
%TED_0EBKH TED_0EBKL 
%TED_0PBKH TED_0PBKL TED_3PBKH TED_3PBKL;

figure(10);
suptitle('SEM: Total Energy Detector');
subplot(3,2,1);
plot(1:numel(TED_0EFL),TED_0EFL,1:numel(TED_0PFL),TED_0PFL,1:numel(TED_3EFL),TED_3EFL,1:numel(TED_3PFL),TED_3PFL);
title('0.05-1 keV Partial Energy Flux');
h_legend=legend('0 Deg Electrons','0 Deg Protons','30 Deg Electrons','30 Deg Protons');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(TED_0EFL) -0.5 256]);

subplot(3,2,2);
plot(1:numel(TED_0EFH),TED_0EFH,1:numel(TED_0PFH),TED_0PFH,1:numel(TED_3EFH),TED_3EFH,1:numel(TED_3PFH),TED_3PFH);
title('2-10 keV Partial Energy Flux');
h_legend=legend('0 Deg Electrons','0 Deg Protons','30 Deg Electrons','30 Deg Protons');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(TED_0EFL) -0.5 256]);

subplot(3,2,3);
plot(1:numel(TED_0DEM),TED_0DEM,1:numel(TED_0DPM),TED_0DPM,1:numel(TED_3DEM),TED_3DEM,1:numel(TED_3DPM),TED_3DPM);
title('Maximum Differential Energy Flux');
h_legend=legend('0 Deg Electrons','0 Deg Protons','30 Deg Electrons','30 Deg Protons');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(TED_0EFL) -0.5 256]);

subplot(3,2,4);
plot(1:numel(TED_0EM),TED_0EM,'-x',1:numel(TED_0PM),TED_0PM,'-x',1:numel(TED_3EM),TED_3EM,'-x',1:numel(TED_3PM),TED_3PM,'-x');
title('Energy of Maximum Differential Energy Flux');
h_legend=legend('0 Deg Electrons','0 Deg Protons','30 Deg Electrons','30 Deg Protons');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(TED_0EFL) -0.5 17]);

subplot(3,2,5);
plot(1:numel(TED_0DP1),TED_0DP1,1:numel(TED_0DP2),TED_0DP2,1:numel(TED_0DP3),TED_0DP3,1:numel(TED_0DP4),TED_0DP4,1:numel(TED_3DP1),TED_3DP1,1:numel(TED_3DP2),TED_3DP2,1:numel(TED_3DP3),TED_3DP3,1:numel(TED_3DP4),TED_3DP4,1:numel(TED_0DP1),TED_0DP1,1:numel(TED_0DP2),TED_0DP2,1:numel(TED_0DP3),TED_0DP3,1:numel(TED_0DP4),TED_0DP4,1:numel(TED_3DP1),TED_3DP1,1:numel(TED_3DP2),TED_3DP2,1:numel(TED_3DP3),TED_3DP3,1:numel(TED_3DP4),TED_3DP4);
%plot(1:numel(TED_0DE1),TED_0DE1,1:numel(TED_0DE2),TED_0DE2,1:numel(TED_0DE3),TED_0DE3,1:numel(TED_0DE4),TED_0DE4,1:numel(TED_3DE1),TED_3DE1,1:numel(TED_3DE2),TED_3DE2,1:numel(TED_3DE3),TED_3DE3,1:numel(TED_3DE4),TED_3DE4);
%plot(1:numel(TED_0DP1),TED_0DP1,1:numel(TED_0DP2),TED_0DP2,1:numel(TED_0DP3),TED_0DP3,1:numel(TED_0DP4),TED_0DP4,1:numel(TED_3DP1),TED_3DP1,1:numel(TED_3DP2),TED_3DP2,1:numel(TED_3DP3),TED_3DP3,1:numel(TED_3DP4),TED_3DP4);
h_legend=legend('0 Deg P1','0 Deg P2','0 Deg P3','0 Deg P4','30 Deg P1','30 Deg P2','30 Deg P3','30 Deg P4','0 Deg E1','0 Deg E2','0 Deg E3','0 Deg E4','30 Deg E1','30 Deg E2','30 Deg E3','30 Deg E4');
set(h_legend,'FontSize',6);
set(h_legend,'Location','best');
title('Four Point Energy/Flux Spectrum');
axis([1 numel(TED_0DP1) -0.5 256]);

subplot(3,2,6);
plot(1:numel(TED_0EBKL),TED_0EBKL,1:numel(TED_0EBKH),TED_0EBKH,1:numel(TED_0PBKL),TED_0PBKL,1:numel(TED_0PBKH),TED_0PBKH);
title('Background');
h_legend=legend('Electron Low','Electron High','Protons Low','Protons High');
set(h_legend,'FontSize',8);
set(h_legend,'Location','best');
axis([1 numel(TED_0EBKL) -0.5 256]);

%% SBUV/2 - Pull Out "Solar Backscatter Ultraviolet Spectral Radiometer" experiment Words
% **IMPORTANT**Usable data can only be collected by the SBUV when it is 
% integrated onto an afternoon spacecraft due to solar angle requirements.

% 3.8.1.1 Purpose of the SBUV Instrument
% The purpose of the SBUV instrument is to measure the solar irradiance and Earth radiance in the near ultraviolet spectrum. From these data, the following atmospheric properties can be deduced:
% 1. The global and vertical distribution of stratospheric ozone
% 2. The structure and dynamics of stratospheric ozone
% 3. Photochemical processes and the influence of "trace" constituents on the ozone layer.
% 4. Long-term solar activity in the Ultraviolet spectrum.


% Bytes 36, 37, 80, 81

%The multiple levels of subcommutation may be justify using ONLY error free
%minor frames for the input data...

clear SBUVLines SBUVdata

SBUVLines(10,1)=0;

minorFrameID = bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6));
SBUVdata(:,1) = bitshift(minorFrames(:,36+1),8)+(minorFrames(:,37+1)); 
SBUVdata(:,2) = bitshift(minorFrames(:,80+1),8)+(minorFrames(:,81+1));

% obsolete data structures below
%SBUVdata(:,3) = (minorFrames(:,80+1));
%SBUVdata(:,4) = (minorFrames(:,81+1));
%Byte 5 is a pseudobyte that holds the minor frame ID corresponding to each
%entry in the previous four items. 
%SBUVdata(:,5) = bitor(bitshift(bitand(minorFrames(:,5),1),8),minorFrames(:,6));

idx = ones(10,1);
for frame=1:numel(minorFrameID)    
    if(mod(minorFrameID(frame)-0,10) == 0) %Minor frames 0, 10, 20, etc
            SBUVLines(idx(1),1,1) = SBUVdata(frame,1); %Status word 1
            SBUVLines(idx(1),1,2) = SBUVdata(frame,2); %Range-1 Data
            SBUVLines(idx(1),1,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(1),1,4) = frame; %MinorFrameID
            idx(1)=idx(1)+1;
    elseif(mod(minorFrameID(frame)-1,10) == 0) %Minor frames 1, 11, 21, etc
            SBUVLines(idx(2),2,1) = SBUVdata(frame,1); %Status word 2
            SBUVLines(idx(2),2,2) = SBUVdata(frame,2); %Range-2 Data
            SBUVLines(idx(2),2,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(2),2,4) = frame; %MinorFrameID
            idx(2)=idx(2)+1;
    elseif(mod(minorFrameID(frame)-2,10) == 0) %Minor frames 2, 12, 22, etc
            SBUVLines(idx(3),3,1) = SBUVdata(frame,1); %Analog Sub Mux
            SBUVLines(idx(3),3,2) = SBUVdata(frame,2); %Range-3 Data
            SBUVLines(idx(3),3,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(3),3,4) = frame; %MinorFrameID
            idx(3)=idx(3)+1;
    elseif(mod(minorFrameID(frame)-3,10) == 0) %Minor frames 3, 13, 23, etc
            SBUVLines(idx(4),4,1) = SBUVdata(frame,1); %Memory Verify
            SBUVLines(idx(4),4,2) = SBUVdata(frame,2); %0x0000
            SBUVLines(idx(4),4,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(4),4,4) = frame; %MinorFrameID
            idx(4)=idx(4)+1;
    elseif(mod(minorFrameID(frame)-4,10) == 0) %Minor frames 4, 14, 24, etc
            SBUVLines(idx(5),5,1) = SBUVdata(frame,1); %Status word 3
            SBUVLines(idx(5),5,2) = SBUVdata(frame,2); %0x0000
            SBUVLines(idx(5),5,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(5),5,4) = frame; %MinorFrameID
            idx(5)=idx(5)+1;
    elseif(mod(minorFrameID(frame)-5,10) == 0) %Minor frames 5, 15, 25, etc
            SBUVLines(idx(6),6,1) = SBUVdata(frame,1); %Status word 4
            SBUVLines(idx(6),6,2) = SBUVdata(frame,2); %0x0000
            SBUVLines(idx(6),6,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(6),6,4) = frame; %MinorFrameID
            idx(6)=idx(6)+1;
    elseif(mod(minorFrameID(frame)-6,10) == 0) %Minor frames 6, 16, 26, etc
            SBUVLines(idx(7),7,1) = SBUVdata(frame,1); %Grating Pos
            SBUVLines(idx(7),7,2) = SBUVdata(frame,2); %0x0000
            SBUVLines(idx(7),7,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(7),7,4) = frame; %MinorFrameID
            idx(7)=idx(7)+1;
    elseif(mod(minorFrameID(frame)-7,10) == 0) %Minor frames 7, 17, 27, etc
            SBUVLines(idx(8),8,1) = SBUVdata(frame,1); %Cloud Cover Radiometer Data
            SBUVLines(idx(8),8,2) = SBUVdata(frame,2); %0x0000
            SBUVLines(idx(8),8,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(8),8,4) = frame; %MinorFrameID
            idx(8)=idx(8)+1;
    elseif(mod(minorFrameID(frame)-8,10) == 0) %Minor frames 8, 18, 28, etc
            SBUVLines(idx(9),9,1) = SBUVdata(frame,1); %Radiomatric DC Level / Grating Pos Error
            SBUVLines(idx(9),9,2) = SBUVdata(frame,2); %0x0000
            SBUVLines(idx(9),9,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(9),9,4) = frame; %MinorFrameID
            idx(9)=idx(9)+1;
    elseif(mod(minorFrameID(frame)-9,10) == 0) %Minor frames 9, 19, 29, etc
            SBUVLines(idx(10),10,1) = SBUVdata(frame,1); %Frame Code Sync
            SBUVLines(idx(10),10,2) = SBUVdata(frame,2); %0x0000
            SBUVLines(idx(10),10,3) = frameTime(frame,36); %Timestamp (of first SBUV byte)
            SBUVLines(idx(10),10,4) = frame; %MinorFrameID
            idx(10)=idx(10)+1;
    end
    %elseif 
end

%% SBUV/2 - Plot Raw Cloud Cover Radiometer Data
%The CCR has a fixed 379 nm filter for wavelength selection and is 
%co-aligned to the monochromator; therefore, it views the same scene as the
%monochromator. The output of the CCR represents the amount of cloud cover 
%in a scene, as the name implies, and is used to remove the effects of 
%clouds in the monochromator data.
stem(SBUVLines(:,7+1,3),SBUVLines(:,7+1,1),'*');
hold on
plot(SBUVLines(:,1+0,3),SBUVLines(:,1+0,1),'-+')
plot(SBUVLines(:,1+5,3),SBUVLines(:,1+5,1),'-x')
hold off

%% SBUV/2 - Cloud Cover Radiometer Dissected Data
% Status word 1 - Measurements issued in pairs with low range first
% followed by high range data
%    bit 11 is the calibration bit which affects the following measurement. 
%    bits 8,9,10 are a counter that count up twice per major frame (unused by this section)

%status word 1 progression, each is issued TWICE in a row
%000000 000 0111010
%000000 001 0111010

%000000 010 0110101
%000000 011 0110101
%000000 100 0110101
%000000 101 0110101
%000000 110 0110101
%000000 111 0110101

%000000 000 0111111
%000000 001 0111111
%000000 010 0111111
%000000 011 0111111

%000001 100 0111111 - Calibration
%000001 101 0111111 - Calibration
%000001 110 0111111 - Calibration
%000001 111 0111111 - Calibration

clear SBUVCCRL SBUVCCRH SBUVCCRCal SBUVCCRLT SBUVCCRHT SBUVCCRCalT;
SBUVCCRL(1) = 0;
SBUVCCRH(1) = 0;
SBUVCCRCal(1) = 0;
SBUVCCRLT(1) = 0;
SBUVCCRHT(1) = 0;
SBUVCCRCalT(1) = 0;

for idx=1:size(SBUVLines,0+1)-1
    if(bitand(SBUVLines(idx,0+1,1),1024) == 1024) %Is cal data?
        %store in cal data place or something
        idx2 = find(SBUVLines(:,7+1,4) > SBUVLines(idx,0+1,4),1); %find the next CCR reading 
        if(idx2 < size(SBUVLines,1))
           SBUVCCRCal(end+1) = SBUVLines(idx2,7+1,1);
           SBUVCCRCalT(end+1) = SBUVLines(idx2,7+1,3);
        end
    elseif(SBUVLines(idx,0+1,1) == SBUVLines(idx+1,0+1,1)) %Is this the first (low range pt) of a pair?
        idx2 = find(SBUVLines(:,7+1,4) > SBUVLines(idx,0+1,4),1); %find the next CCR reading 
        if(idx2 < size(SBUVLines,1))
           SBUVCCRL(end+1) = SBUVLines(idx2,7+1,1);
           SBUVCCRLT(end+1) = SBUVLines(idx2,7+1,3);
        end
        %store low range array
    elseif(SBUVLines(idx,0+1,1) == SBUVLines(idx-1,0+1,1)) %Is this the second (high range pt) of a pair?
        idx2 = find(SBUVLines(:,7+1,4) > SBUVLines(idx,0+1,4),1); %find the next CCR reading 
        if(idx2 < size(SBUVLines,1))
           SBUVCCRH(end+1) = SBUVLines(idx2,7+1,1);
           SBUVCCRHT(end+1) = SBUVLines(idx2,7+1,3);
        end
        %store high range array
    end
end
plot(SBUVCCRHT,SBUVCCRH,'-x',SBUVCCRLT,SBUVCCRL,'-x');
%,SBUVCCRCalT,SBUVCCRCal,'-x');

%find(SBUVLines(:,7+1,4) > 600,1) %find next frame that corresponds to this


%% SBUV/2 - Plot Raw Telemetry Channels and Channel Counter
%plot(SBUVLines(:,2+1,3),SBUVLines(:,2+1,1),'-* ',SBUVLines(:,9+1,3),bitshift(SBUVLines(:,9+1,1),-8),'*');
[hAx,hLine1,hLine2] = plotyy(SBUVLines(:,2+1,3),bitshift(SBUVLines(:,2+1,1),-8),SBUVLines(:,9+1,3),bitshift(SBUVLines(:,9+1,1),-12));
hLine1.Marker = 'x';
hLine2.Marker = '*';

%% SBUV/2 - Plot grating position and some status crap and some other crap too
plot(SBUVLines(:,0+1,3),SBUVLines(:,0+1,2),'-x',SBUVLines(:,1+1,3),SBUVLines(:,1+1,2),'-o',SBUVLines(:,2+1,3),SBUVLines(:,2+1,2),'-+');
hold on
stem(SBUVLines(:,6+1,3),SBUVLines(:,6+1,1),'*');
plot(SBUVLines(:,0+2,3),SBUVLines(:,0+2,1),'-+')
hold off

plot(SBUVLines(:,1,2),'-x') %range data 1
hold on
plot(SBUVLines(:,2,2),'-o') %range data 2
plot(SBUVLines(:,3,2),'-+') %range data 3
plot(SBUVLines(:,6+1,1),'-+') %grating pos
hold off