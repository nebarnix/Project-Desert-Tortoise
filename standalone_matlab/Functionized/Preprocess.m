%% *Load Wav FIle
%tic;
clear all;
%hfile = 'C:\Users\nebarnix\Documents\vmshare\POES\35-N19\15-53-35_137770kHz_50k.wav';
hfile = 'C:\Users\nebarnix\Documents\vmshare\POES\33-N18 hackrf\18-01-50_137350kHz_50k.wav';
%hfile = 'C:\sdrsharp\sdr-install 1347\sdrsharp\IF\2016_04_25\07-10-40_137350kHz.wav'; 
%fid = fopen(hfile,'rb');
[audioData,Fs] = audioread(hfile);
%mmInDataStream = fread(fid,'float32');

%dataStreamIn = (audioData(:,1) + 1i*audioData(:,2))';

%swap I/Q
dataStreamIn = (audioData(:,2) + 1i*audioData(:,1))';

Ts = 1/(Fs);
BasebandRawTime=0:Ts:Ts*(numel(dataStreamIn)-1);
%spectrogram(dataStreamIn(22207500:22257500),blackman(128),60,128,1e3)
n=1;
fprintf(['Loaded IQ WAV file, Sample rate detected as ' num2str(Fs/1000) 'Ksps\n']);

%% *Load Raw File
%tic;
clear all;
%hfile = 'C:\Users\nebarnix\Documents\vmshare\POES\6-NOAA15\baseband\NOAA15_06092015_telemetry_50k.raw'; 
hfile = 'F:\POES\6-NOAA15\NOAA15_06092015_telemetry_50k.raw';
 
fid = fopen(hfile,'rb');
dataStreamIn = fread(fid,'float32');
fclose(fid);

dataStreamIn = (dataStreamIn(2:2:end) + 1i*dataStreamIn(1:2:end))';

Nfft = 1500000;

Fs=Fs;
Ts = 1/(Fs);
BasebandRawTime=0:Ts:Ts*(numel(dataStreamIn)-1);
n=1;
fprintf(['Loaded IQ RAW file, using ' num2str(Fs/1000) 'Ksps as sample rate\n']);
%% Process Data 

% %% Plot spectrogram of input data, piecewise
% Figure(1);
% inc = 1000000;
% n=n+inc; 
% spectrogram(dataStreamIn(n:inc+n),blackman(128),60,32,1e3);

%static gain block, use the middle of the recording, set to sin wave
%average value of 0.637
dataStreamIn = StaticGain(dataStreamIn(end/2:end/2+10000),10000,0.637).*dataStreamIn;

%[dataStreamAGC] = NormalizingAGC(dataStreamIn,.0000001);
%[dataStreamAGC] = NormalizingAGC(dataStreamIn,.00001);

%PLL Carrier Tracking Loop
%clear dataStreamPLLOut;
[dataStreamPLLOut, d_freqi, d_locksigi, firstLock] = CarrierTrackPLL(dataStreamIn, Fs, 4500, 0.1, 0.01, 0.001);
fprintf(['Locked at ' num2str(d_freqi(firstLock)) 'Hz\n']);
% %% Write results to wav file, making sure no value exceeds 1...
% audiowrite('PLL.wav',dataStreamPLLOut ./ max(dataStreamPLLOut),Fs,'BitsPerSample',16);

% %Plot Lock Stuffs piecewise
% figure(2);
% subplot(2,1,1);
% spectrogram(dataStreamPLLOut(100000:200000),blackman(128),60,128,1e3)
% subplot(2,1,2);
% plot(d_locksigi(1:200000));
% axis([0 200000 -.001 .015]);
% n=1;

% %Plot Whole Spectrogram
% figure(3);
% spectrogram(dataStreamPLLOut(1:end),blackman(128),60,128,Fs)
% %axis([0 200000 -.001 .015]);
% n=1;

% %Plot lock signal and spectrogram piecewise
% figure(4);
% n=n+200000; 
% subplot(2,1,1);
% spectrogram(dataStreamPLLOut(n:200000+n),blackman(128),60,128,Fs);
% subplot(2,1,2);
% plot(d_locksigi(n:200000+n));
% axis('tight');
% %axis([0 200000 -.001 .015]);

 %Plot lock and frequency
 figure(5);
 plotyy(1:numel(d_locksigi(1:10000:end)),d_locksigi(1:10000:end),1:numel(d_freqi(1:10000:end)),d_freqi(1:10000:end).*Fs./(2.*pi));
 %axis 'tight';

% %Plot bit strengths
% figure(6);
% %subplot(2,1,1);
% pointScale = 51;
% %length(y(1:pointScale:end));
% %scatter(mmRawTime(1:pointScale:end),real(mmInDataStream(1:pointScale:end)),2,'.');
% %plot(RawTime(1:pointScale:end),imag(y(1:pointScale:end)),2,'.');
% plot(dataStreamPLLOut(1:51:end),'.');
% axis([0 max(PLLRawTime) -2 2]);

%Lowpass filter the data
%clear dataStreamLPF;

fprintf('Lowpass Filtering...');
%remez(42,[0 0.37 0.43 1],[1 1 0 0]);
%dataStreamLPF = filter(LPF,1,dataStreamAGC);
filterOrder = 26;
%LPF = fir1(filterOrder,(11e3/(Fs/2)),kaiser(filterOrder+1,6.76));
filterTaps = fir1(filterOrder,(11e3/(Fs/2)),blackman(filterOrder+1));
%LPF = fir1(filterOrder,(11e3/(Fs/2)));
dataStreamLPF = filter(filterTaps,1,dataStreamPLLOut);
%dataStreamLPF = filter(1/4*[1 1 1 1],1,dataStreamAGC); %moving avg
fprintf('done.\n');

% %Plot 
% figure(7);
% subplot(1,2,1);
% spectrogram(dataStreamPLLOut(1:1e5),blackman(128),60,128,Fs);
% subplot(1,2,2);
% spectrogram(dataStreamLPF(1:1e5),blackman(128),60,128,Fs);

%Automatic Gain Control Block
%clear dataStreamAGC;
[dataStreamAGC, gaini] = NormalizingAGC(dataStreamLPF,.00025);

% %Save AGC wav file, normalize max value to 1
% audiowrite('PLL_AGC.wav',dataStreamAGC/max(dataStreamAGC),Fs,'BitsPerSample',16);

% %Plot
% figure(8);
% subplot(2,1,1)
% plot(PLLRawTime,dataStreamLPF,PLLRawTime,dataStreamAGC);
% subplot(2,1,2)
% plot(PLLRawTime,d_locksigi,PLLRawTime,gaini);
% %[hAx,hLine1,hLine2] = plotyy(RawTime(1:num_loop),dataStreamLPF(1:num_loop),RawTime(1:num_loop),dataStreamAGC(1:num_loop));
% %hLine1.LineStyle = '--';
% %hLine2.LineStyle = ':';
% %axis([1 10000 -1.5 1.5 ]);

% *MM Clock Recovery*
%[dataStreamOut, dataStreamOutTime] = UpsamplingMMClockRecovery(dataStreamIn, FsIn, FsOut, baud, stepSpread);
%clear dataStreamBits dataStreamBitsTime;
[dataStreamBits, dataStreamBitsTime] = UpsamplingMMClockRecovery(dataStreamAGC, BasebandRawTime, Fs, 8320*15, 8320*2-1, 10);
%[dataStreamOut, dataStreamOutTime] = MMClockRecovery(dataStreamIn, Fs, baud, stepSpread);
%[dataStreamBits, dataStreamBitsTime] = MMClockRecovery(dataStreamAGC, Fs, 8320*2, 10);

% %Plot clocksync'd bit strengths
% figure(9);
% %subplot(2,1,2);
% pointScale = 1;
% 
% %startT = 365;
% %endT = 375;
% 
% startT = 320;
% endT = 325;
% 
% startRawT=find(mmRawTime >= startT,1);
% startBitT=find(mmbitTime >= startT,1);
% endRawT=find(mmRawTime >= endT,1);
% endBitT=find(mmbitTime >= endT,1);
% %length(y(1:pointScale:end));
% plot(mmRawTime(startRawT:pointScale:endRawT),real(mmInDataStream(startRawT:pointScale:endRawT)),'.-',mmbitTime(startBitT:pointScale:endBitT),real(mmOutBits(startBitT:pointScale:endBitT)),'x');
% %hold on;
% %scatter(bitTime(1:pointScale:end),real(Bits(1:pointScale:end)),2,'.');
% 
% %hold off;
% %plot(RawTime(1:pointScale:end),imag(y(1:pointScale:end)),2,'.');
% axis([mmbitTime(startBitT) mmbitTime(endBitT) -10 10]);
% toc;

%Run functionized

%manchester threshold to the inverse gain of what it takes to get to 1 which should make
%it right at the center of the constellation points
[machesterStreamBits, machesterStreamBitime] = manchesterDecodeFloat(dataStreamBits, dataStreamBitsTime, 1/StaticGain(dataStreamBits(end/2:end/2+10000),10000,1));
[SyncWordIndex, SyncWordInvIndex] = syncWordDetect(machesterStreamBits);
[minorFrames, frameTime] = convertBitsToBytes(machesterStreamBits, machesterStreamBitime, SyncWordIndex, SyncWordInvIndex);
[dayNum, spaceCraft, minorFrameID] = daytimeDecode(minorFrames, frameTime);
%plotyy(BasebandRawTime(1:1000:end),d_locksigi(1:1000:end),frameTime(:,5),minorFrameID);
[goodFrames, parity] = checkParity(minorFrames);

%% Process Data no fuss (no generated plots)

dataStreamIn = StaticGain(dataStreamIn(end/2:end/2+10000),10000,0.637).*dataStreamIn;

[dataStreamPLLOut, d_freqi, d_locksigi, firstLock] = CarrierTrackPLL(dataStreamIn, Fs, 4500, 0.025, 0.0075*2, 0.005/4); %Optimized at 50ksps, probably changes at different sample rates?
fprintf(['Locked at ' num2str(d_freqi(firstLock)*Fs/(2*pi)) 'Hz\n']);

%Plot lock and frequency
figure(10);
plotyy(1:numel(d_locksigi(1:10000:end)),d_locksigi(1:10000:end),1:numel(d_freqi(1:10000:end)),d_freqi(1:10000:end).*Fs./(2.*pi));
%Lowpass filter the data

fprintf('Lowpass Filtering...');
%remez(42,[0 0.37 0.43 1],[1 1 0 0]);
%dataStreamLPF = filter(LPF,1,dataStreamAGC);
filterOrder = 26; %26 was the default previously
%LPFTaps = fir1(filterOrder,(11e3/(Fs/2)),kaiser(filterOrder+1,6.76));
%LPFTaps = fir1(filterOrder,(11e3/(Fs/2)),blackman(filterOrder+1));
filterTaps = fir1(filterOrder,(11e3/(Fs/2)),blackman(filterOrder+1));
%LPFTaps = fir1(filterOrder,(11e3/(Fs/2)));
dataStreamLPF = filter(filterTaps,1,dataStreamPLLOut);
%dataStreamLPF = filter(1/4*[1 1 1 1],1,dataStreamAGC); %moving avg
fprintf('done.\n');

[dataStreamAGC, gaini] = NormalizingAGC(dataStreamLPF,.00025);


[dataStreamBits, dataStreamBitsTime] = UpsamplingMMClockRecovery(dataStreamAGC, BasebandRawTime, Fs, 8320*18, 8320*2-1, 10, 0.15);

fprintf('Decoding Manchester Bits...'); [machesterStreamBits, machesterStreamBitime] = manchesterDecodeFloat(dataStreamBits, dataStreamBitsTime, 1/StaticGain(dataStreamBits(end/2:end/2+10000),10000,1)); %fprintf('Manchester Decoding Done.\n');
fprintf('Detecting Syncword...'); [SyncWordIndex, SyncWordInvIndex] = syncWordDetect(machesterStreamBits); %fprintf('Syncword Detect Done.\n');
fprintf('Converting Bits to Frames...'); [minorFrames, frameTime] = convertBitsToBytes(machesterStreamBits, machesterStreamBitime, SyncWordIndex, SyncWordInvIndex); fprintf('done.\n');
fprintf('Decoding Day and Time...\n'); [dayNum, spaceCraft, minorFrameID] = daytimeDecode(minorFrames, frameTime);  %fprintf('Day Time Decoding Finished.\n');
%plotyy(BasebandRawTime(1:1000:end),d_locksigi(1:1000:end),frameTime(:,5),minorFrameID);
fprintf('Checking Parity...'); [goodFrames, parity] = checkParity(minorFrames);  %fprintf('Parity Check Done.\n');
yield(sampleInterpt,1) = goodFrames;
yield(sampleInterpt,2) = size(minorFrames,1);

%plot(1:size(yield,1),yield(:,1),1:size(yield,1),yield(:,2))
