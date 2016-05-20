#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "CarrierTrackPLL.h"

double CarrierTrackPLL(double complex *complexDataIn, double *realDataOut, unsigned int nSamples, double Fs, double freqRange, double d_lock_threshold, double loopbw_acq, double loopbw_track)
   {
   double bw = loopbw_acq;
   double sample_phase;
   
   static double firstLock = -1;
   
   static double damp=1;
   static double d_alpha;
   static double d_beta;
      
   static double d_phase;
   static double d_freq;
   static double d_max_freq;
   static double d_min_freq;
   
   if(firstLock == -1)
      {     
      d_alpha = (4 * damp * bw) / (1 + 2 * damp * bw + bw * bw);
      d_beta = (4 * bw * bw) / (1 + 2 * damp * bw + bw * bw);
      
      d_phase = 0.1; //something not zero for benchmarking the function (sin(0) is probably a shortcut)
      d_freq  = 2.0*M_PI*0 / Fs;
      d_max_freq   = 2.0*M_PI*freqRange/Fs; //+/-4500 for 2m polar sats
      d_min_freq   = -2.0*M_PI*freqRange/Fs;
      firstLock = 0;
      }
   static double d_locksig = 0;
   static double lockSigAlpha = 0.00005;
   
   double t_imag;
   double t_real;
   double error;
   unsigned int idx;
   //dataStreamOut = zeros(1,size(dataStreamIn,2));
   //d_freqi = zeros(1,size(dataStreamIn,2));
   //d_locksigi = zeros(1,size(dataStreamIn,2));
   
   //double progress = 0;
   //double percentcomplete=0;
   //double onePercent = numel(dataStreamOut) / 100;
   
   //printf("Carrier Tracking PLL:");
   //MSG = [convertnum(percentcomplete) '%'];
   //msgCount = numel(MSG);
   //printf("%s",MSG);
   
           
   for (idx=0; idx<nSamples; idx++)
      {
       //sincos(d_phase, &t_imag, &t_real);
       
       /*progress = progress+1;
       if(progress >= 1*onePercent)        
           for progress=1:msgCount
              fprintf('%c',char(8));
           end
           percentcomplete = percentcomplete + 1;
           MSG = [convertnum(percentcomplete) '%'];
           //MSG = sprintf('%d%%',percentcomplete);
           msgCount = numel(MSG);
           fprintf('%s',MSG);        
           progress = 0;
       end*/
       
       t_imag = sin(d_phase);
       t_real = cos(d_phase);    
       
       //t_imag = sin_lookup(d_phase);
       //t_real = cos_lookup(d_phase);
       
       //t_imag = sinlutydata(floor((d_phase - (-2*pi))/((2*pi - -2*pi)/1024)));    
       //t_real = coslutydata(floor((d_phase - (-2*pi))/((2*pi - -2*pi)/1024)));    
       
       //t_imag = sinlutydata(floor((d_phase +2*pi)/((12.5664)/nlut)));    
       //t_real = coslutydata(floor((d_phase +2*pi)/((12.5664)/nlut)));    
       
   
       //t_imag = coslutydata(coslutxdata(floor(d_phase*1024)));
       //t_real = sinlutydata(sinlutxdata(floor(d_phase*1024)));
       
       //Shift frequency by loop phase
       
       realDataOut[idx] = cimagf(complexDataIn[idx] * (t_real+I*-t_imag) );
       
       //double re, im;
       //gr::sincosf(d_phase, &im, &re);
       //out[i] = (in[i]*gr_complex(re, -im)).imag();
       
       //Calculate Error
       //sample_phase = atan2_approx(imag(dataStreamIn(idx)),real(dataStreamIn(idx)));
       sample_phase = atan2(cimag(complexDataIn[idx]),creal(complexDataIn[idx]));
       
       //error = mod_2pi(sample_phase-d_phase);
       
       //mod2pi the error
       if((sample_phase-d_phase) > M_PI)
           error = (sample_phase-d_phase)-2*M_PI;
       else if((sample_phase-d_phase) < -M_PI)
           error = (sample_phase-d_phase)+2*M_PI;
       else
           error= sample_phase-d_phase;
       
       
       //advance the loop
       d_freq = d_freq + d_beta * error;
       d_phase = d_phase + d_freq + d_alpha * error;
       
       //wrap the phase
       while(d_phase > 2*M_PI)
           d_phase = d_phase-2*M_PI;
       
       while(d_phase < -2*M_PI)
           d_phase = d_phase+2*M_PI;
       
       //Limit the frequency
       if(d_freq > d_max_freq)
           d_freq = d_max_freq;
       else if(d_freq < d_min_freq)
           d_freq = d_min_freq;
       
       
      ///////////d_freqi(idx) = d_freq;
      
      //d_locksig = d_locksig * (1.0 - d_alpha) + d_alpha*(real(dataStreamIn(idx)) * t_real + imag(dataStreamIn(idx)) * t_imag);
      d_locksig = d_locksig * (1.0 - lockSigAlpha) + lockSigAlpha*(creal(complexDataIn[idx]) * t_real + cimag(complexDataIn[idx]) * t_imag);
       
      //moving average filter the locksig. For loops are slow in matlab.
      //Circshift instead?
      
      //It isn't actually neccesary to preserve order with a rectangular window. 
      //You can just roll the location of the placed value 1-10?
      //for idx2=0:locksigMovAvgOrder-2        
      //   //fprintf('%d to //d\n',locksigMovAvgOrder-idx2,locksigMovAvgOrder-idx2-1);
      //   locksigMovAvg(locksigMovAvgOrder-idx2) = locksigMovAvg(locksigMovAvgOrder-idx2-1);
      //end
      //locksigMovAvg = circshift(locksigMovAvg,1,2); //shift the array once down the line
      
      //locksigMovAvg(mod(idx,locksigMovAvgOrder)+1) = d_locksig; 
      //d_locksigi(idx) = sum(locksigMovAvg)/locksigMovAvgOrder; //average
      //lockSignalAccumulator = (lockSigAlpha * d_locksig) + (1.0 - lockSigAlpha) * lockSignalAccumulator;
      
      ///////////////d_locksigi(idx) = d_locksig;
      
      //d_locksigi(idx) = mean(locksigMovAvg); //faster? SLOWER! eep!
      
      //d_locksig > d_lock_threshold (0.01)
      if(d_locksig > d_lock_threshold && firstLock == 0) //This needs to be a moving average or at least somewhat smoothed
         {
         printf(" : PLL locked at %0.2fHz\n", d_freq*Fs/(2.0*M_PI));
         //fprintf(['PLL locked at ' num2str(d_freq*Fs/(2*pi)) '\n']);
         firstLock = idx;
         //d_alpha = 0.001; //decrease to tracking gain
         bw = loopbw_track;        
         d_alpha = (4 * damp * bw) / (1 + 2 * damp * bw + bw * bw);
         d_beta = (4 * bw * bw) / (1 + 2 * damp * bw + bw * bw);
         //d_alpha = alpha_track;
         //d_beta = d_alpha/25;
         }
      }
   //printf("%f\t%f\n",d_locksig,d_freq*Fs/(2.0*M_PI));   
   return firstLock;
   }
/*
for progress=1:msgCount
    fprintf('%c',char(8));        
end
fprintf('100%%\n');

}*/
/*
function s = convertnum(n)
   s = [];
   while n > 0
      d = mod(n,10);
      s = [char(48+d), s];
      n = (n-d)/10;
   end
end*/