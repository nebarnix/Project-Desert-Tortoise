%#codegen
%PLL Carrier Tracking Loop
function [dataStreamOut, d_freqi, d_locksigi, firstLock] = CarrierTrackPLL(dataStreamIn, Fs, freqRange, d_lock_threshold, loopbw_acq, loopbw_track)
%d_lock_threshold = 0.5;
%d_alpha = 0.01; %aquisition gain
%d_alpha = loopbw_acq;
%d_beta = d_alpha/25;
bw = loopbw_acq;
damp=1;
d_alpha = (4 * damp * bw) / (1 + 2 * damp * bw + bw * bw);
d_beta = (4 * bw * bw) / (1 + 2 * damp * bw + bw * bw);

d_phase = 0.1; %something not zero for benchmarking the function (sin(0) is probably a shortcut)
d_freq  = 2*pi*1267/Fs;
d_max_freq   = 2*pi*freqRange/Fs; %+/-4500 for 2m polar sats
d_min_freq   = -2*pi*freqRange/Fs;
d_locksig = 0;
firstLock = 0;
dataStreamOut = zeros(1,size(dataStreamIn,2));
d_freqi = zeros(1,size(dataStreamIn,2));
d_locksigi = zeros(1,size(dataStreamIn,2));
%nlut = 1024;
%[coslutxdata, coslutydata] = create_cos_lookup(nlut);
%[sinlutxdata, sinlutydata] = create_sin_lookup(nlut);

progress = 0;
percentcomplete=0;
onePercent = numel(dataStreamOut) / 100;
%msgCount = 1;
fprintf('Carrier Tracking PLL:');
MSG = [convertnum(percentcomplete) '%'];
msgCount = numel(MSG);
fprintf('%s',MSG);

lockSigAlpha = 0.00005;
        
for idx=1:numel(dataStreamIn)
    %sincosf(d_phase, &t_imag, &t_real);
    
    progress = progress+1;
    if(progress >= 1*onePercent)        
        for progress=1:msgCount
           fprintf('%c',char(8));
        end
        percentcomplete = percentcomplete + 1;
        MSG = [convertnum(percentcomplete) '%'];
        %MSG = sprintf('%d%%',percentcomplete);
        msgCount = numel(MSG);
        fprintf('%s',MSG);        
        progress = 0;
    end
    
    t_imag = sin(d_phase);
    t_real = cos(d_phase);    
    
    %t_imag = sin_lookup(d_phase);
    %t_real = cos_lookup(d_phase);
    
    %t_imag = sinlutydata(floor((d_phase - (-2*pi))/((2*pi - -2*pi)/1024)));    
    %t_real = coslutydata(floor((d_phase - (-2*pi))/((2*pi - -2*pi)/1024)));    
    
    %t_imag = sinlutydata(floor((d_phase +2*pi)/((12.5664)/nlut)));    
    %t_real = coslutydata(floor((d_phase +2*pi)/((12.5664)/nlut)));    
    

    %t_imag = coslutydata(coslutxdata(floor(d_phase*1024)));
    %t_real = sinlutydata(sinlutxdata(floor(d_phase*1024)));
    
    %Shift frequency by loop phase
    
    dataStreamOut(idx) = imag(dataStreamIn(idx) * (t_real+1i*-t_imag));
    
    %float re, im;
    %gr::sincosf(d_phase, &im, &re);
    %out[i] = (in[i]*gr_complex(re, -im)).imag();
    
    
    %Calculate Error
    %sample_phase = atan2_approx(imag(dataStreamIn(idx)),real(dataStreamIn(idx)));
    sample_phase = atan2(imag(dataStreamIn(idx)),real(dataStreamIn(idx)));
    
    %error = mod_2pi(sample_phase-d_phase);
    
    %mod2pi the error
    if((sample_phase-d_phase) > pi)
        error= (sample_phase-d_phase)-2*pi;
    elseif((sample_phase-d_phase) < -pi)
        error = (sample_phase-d_phase)+2*pi;
    else
        error= sample_phase-d_phase;
    end
    
    %advance the loop
    d_freq = d_freq + d_beta * error;
    d_phase = d_phase + d_freq + d_alpha * error;
    
    %wrap the phase
    while(d_phase > 2*pi)
        d_phase = d_phase-2*pi;
    end
    
    while(d_phase < -2*pi)
        d_phase = d_phase+2*pi;
    end
    
    %Limit the frequency
    if(d_freq > d_max_freq)
        d_freq = d_max_freq;
    elseif(d_freq < d_min_freq)
        d_freq = d_min_freq;
    end
    
    d_freqi(idx) = d_freq;
    
    %d_locksig = d_locksig * (1.0 - d_alpha) + d_alpha*(real(dataStreamIn(idx)) * t_real + imag(dataStreamIn(idx)) * t_imag);
    d_locksig = d_locksig * (1.0 - lockSigAlpha) + lockSigAlpha*(real(dataStreamIn(idx)) * t_real + imag(dataStreamIn(idx)) * t_imag);
    
    %moving average filter the locksig. For loops are slow in matlab.
    %Circshift instead?
    
    %It isn't actually neccesary to preserve order with a rectangular window. 
    %You can just roll the location of the placed value 1-10?
    %for idx2=0:locksigMovAvgOrder-2        
    %   %fprintf('%d to %d\n',locksigMovAvgOrder-idx2,locksigMovAvgOrder-idx2-1);
    %   locksigMovAvg(locksigMovAvgOrder-idx2) = locksigMovAvg(locksigMovAvgOrder-idx2-1);
    %end
    %locksigMovAvg = circshift(locksigMovAvg,1,2); %shift the array once down the line
    
    %locksigMovAvg(mod(idx,locksigMovAvgOrder)+1) = d_locksig; 
    %d_locksigi(idx) = sum(locksigMovAvg)/locksigMovAvgOrder; %average
    %lockSignalAccumulator = (lockSigAlpha * d_locksig) + (1.0 - lockSigAlpha) * lockSignalAccumulator;
    d_locksigi(idx) = d_locksig;
    %d_locksigi(idx) = mean(locksigMovAvg); %faster? SLOWER! eep!
    
    
    
    
    %d_locksig > d_lock_threshold (0.01)
    if(d_locksig > d_lock_threshold && firstLock == 0) %This needs to be a moving average or at least somewhat smoothed
        %fprintf(['PLL locked at ' num2str(d_freq*Fs/(2*pi)) '\n']);
        firstLock = idx;
        %d_alpha = 0.001; %decrease to tracking gain
        bw = loopbw_track;        
        d_alpha = (4 * damp * bw) / (1 + 2 * damp * bw + bw * bw);
        d_beta = (4 * bw * bw) / (1 + 2 * damp * bw + bw * bw);
        %d_alpha = alpha_track;
        %d_beta = d_alpha/25;
    end
end
for progress=1:msgCount
    fprintf('%c',char(8));        
end
fprintf('100%%\n');
end

function s = convertnum(n)
   s = [];
   while n > 0
      d = mod(n,10);
      s = [char(48+d), s];
      n = (n-d)/10;
   end
end