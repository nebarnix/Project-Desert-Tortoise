#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "CarrierTrackPLL.h"

double CarrierTrackPLL(double complex *complexDataIn, double *realDataOut, double *lockSignalStreamOut, unsigned int nSamples, double Fs, double freqRange, double d_lock_threshold, double lockSigAlpha, double loopbw_acq, double loopbw_track)
   {
   double bw = loopbw_acq;
   double sample_phase;
   
   static long firstLock = -2;
   
   static double damp=0.999;
   static double d_alpha;
   static double d_beta;
      
   static double d_phase;
   static double d_freq;
   static double d_max_freq;
   static double d_min_freq;
   
   static double averagePhase; 
   
   static double d_locksig = 0;
   //static double lockSigAlpha = 0.00005;
   
   double t_imag;
   double t_real;
   double PLLOutSamplePhase,error;
   double averagePhaseAlpha=0.00005;
   unsigned int idx;
   
   complex double PLLOutSample;
   
   if(firstLock == -2)
      {     
      d_alpha = (4 * damp * bw) / (1 + 2 * damp * bw + bw * bw);
      d_beta = (4 * bw * bw) / (1 + 2 * damp * bw + bw * bw);
      
      d_phase = 0.1; //something not zero for benchmarking the function (sin(0) is probably a shortcut)
      d_freq  = 2.0*M_PI*0 / Fs;
      d_max_freq   = 2.0*M_PI*freqRange/Fs; //+/-4500 for 2m polar sats
      d_min_freq   = -2.0*M_PI*freqRange/Fs;
      firstLock = -1;
      averagePhase=1.5708;
      }
           
   for (idx=0; idx<nSamples; idx++)
      {
       
       t_imag = sin(d_phase);
       t_real = cos(d_phase);    
       
       //shift the frequency by the carrier 
       PLLOutSample = complexDataIn[idx] * (t_real+I*-t_imag);
       
       //data bits are in the imaginary part
       realDataOut[idx] = cimag(PLLOutSample);
       
       
       //calculate phase angle for quality estimation 
       PLLOutSamplePhase = atan2(cimag(PLLOutSample),creal(PLLOutSample));
       
       //running average of the absolute value of the phase
       averagePhase = averagePhase * (1.0 - averagePhaseAlpha) + averagePhaseAlpha * fabs(PLLOutSamplePhase);
       
       //Calculate Error
       sample_phase = atan2(cimag(complexDataIn[idx]),creal(complexDataIn[idx]));
       
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
       
      //d_locksig = d_locksig * (1.0 - d_alpha) + d_alpha*(real(dataStreamIn(idx)) * t_real + imag(dataStreamIn(idx)) * t_imag);
      d_locksig = d_locksig * (1.0 - lockSigAlpha) + lockSigAlpha*(creal(complexDataIn[idx]) * t_real + cimag(complexDataIn[idx]) * t_imag);
      
      if(lockSignalStreamOut != NULL)  
          lockSignalStreamOut[idx] = d_locksig;
      
      //d_locksig > d_lock_threshold (0.01)
      if(d_locksig > d_lock_threshold && firstLock == -1)
         {
         printf(" : PLL locked at %0.2fHz\n", d_freq*Fs/(2.0*M_PI));
         firstLock = idx;
         bw = loopbw_track;        
         d_alpha = (4.0 * damp * bw) / (1.0 + 2.0 * damp * bw + bw * bw);
         d_beta = (4.0 * bw * bw) / (1.0 + 2.0 * damp * bw + bw * bw);
         }
      }
   //printf("%f\t%f\n",d_locksig,d_freq*Fs/(2.0*M_PI));   
   return averagePhase;
   }