#include <complex.h>
#include <tgmath.h>
#include <stdio.h>
#include "CarrierTrackPLL.h"
//#include "atan2_approx.h"

#define RECIP_ITER 1
#define USE_ATAN_APPROX 1

#define coeff_1 (0.78539816339744825)
#define coeff_2 (2.35619449019234475)

//-----------------------------------------------
// Fast arctan2
DECIMAL_TYPE arctan2(DECIMAL_TYPE y, DECIMAL_TYPE x)
{
   DECIMAL_TYPE r, angle, abs_y;

   
   #if USE_FLOATS==1
      abs_y = fabsf(y)+1e-10;      // kludge to prevent 0/0 condition
   #else 
      abs_y = fabs(y)+1e-10;      // kludge to prevent 0/0 condition
   #endif  
   
   if (x >= 0)
      {
      r = (x - abs_y) / (x + abs_y);
      angle = coeff_1 - coeff_1 * r;
      }
   else
      {
      r = (x + abs_y) / (abs_y - x);
      angle = coeff_2 - coeff_1 * r;
      }
   if (y < 0)
      return(-angle);     // negate if in quad III or IV
   else
      return(angle);
}

//could this also work for a double? (with modification maybe)
float Q_rsqrt( float x )
{
float xhalf = 0.5f * x;
int i = *(int*)&x;            // store floating-point bits in integer
i = 0x5f3759df - (i >> 1);    // initial guess for Newton's method
x = *(float*)&i;              // convert new bits into float
x = x*(1.5f - xhalf*x*x);     // One round of Newton's method
x = x*(1.5f - xhalf*x*x);     // Two rounds of Newton's method
return x;
}

DECIMAL_TYPE CarrierTrackPLL(DECIMAL_TYPE complex *complexDataIn, DECIMAL_TYPE *realDataOut, DECIMAL_TYPE *lockSignalStreamOut, unsigned int nSamples, DECIMAL_TYPE Fs, DECIMAL_TYPE freqRange, DECIMAL_TYPE d_lock_threshold, DECIMAL_TYPE lockSigAlpha, DECIMAL_TYPE loopbw_acq, DECIMAL_TYPE loopbw_track)
   {
   
   DECIMAL_TYPE bw = loopbw_acq;
   DECIMAL_TYPE sample_phase;
   
   static long firstLock = -2;
   
   static DECIMAL_TYPE damp=0.999;
   static DECIMAL_TYPE d_alpha;
   static DECIMAL_TYPE d_beta;
      
   static DECIMAL_TYPE d_phase;
   static DECIMAL_TYPE d_freq;
   static DECIMAL_TYPE d_max_freq;
   static DECIMAL_TYPE d_min_freq;
   
   static DECIMAL_TYPE averagePhase; 
   
   static DECIMAL_TYPE d_locksig = 0;
   //static double lockSigAlpha = 0.00005;
   static DECIMAL_TYPE sweep; //radians per second
   
   DECIMAL_TYPE t_imag;
   DECIMAL_TYPE t_real;
   DECIMAL_TYPE PLLOutSamplePhase,error;
   DECIMAL_TYPE averagePhaseAlpha=0.00005;
   unsigned int idx;
   
   DECIMAL_TYPE complexDataInRe,complexDataInIm,invsqrt,magnitude_squared;
   
   //complex double normalizedDataIn;
   complex DECIMAL_TYPE PLLOutSample;
   
   if(firstLock == -2) //Initialize PLL if first run through
      {     
      d_alpha = (4 * damp * bw) / (1 + 2 * damp * bw + bw * bw);
      d_beta = (4 * bw * bw) / (1 + 2 * damp * bw + bw * bw);
      
      d_phase = 0.1; //something not zero for benchmarking the function (sin(0) is probably a shortcut)
      d_freq  = 2.0*M_PI*0 / Fs; //start at 0Hz because we don't know if IQ or QI (we hope we know, but we really don't)
      d_max_freq   = 2.0*M_PI*freqRange/Fs; //+/-4500 for 2m polar sats
      d_min_freq   = -2.0*M_PI*freqRange/Fs;
      firstLock = -1;
      averagePhase=M_PI/2.0;
      sweep = 0.2*(2.0*M_PI/Fs); //0.2
      }
           
   for (idx=0; idx<nSamples; idx++)
      {
       
      #if USE_FLOATS==1
         t_imag = sinf(d_phase);
         t_real = cosf(d_phase);
         
         //shift the frequency by the carrier 
         PLLOutSample = complexDataIn[idx] * (t_real+I*-t_imag);
          
         //data bits are in the imaginary part
         realDataOut[idx] = cimagf(PLLOutSample);
          
         //calculate phase angle for quality estimation 
         #if USE_ATAN_APPROX==1
         PLLOutSamplePhase = arctan2(cimagf(PLLOutSample),crealf(PLLOutSample));
         #else
         PLLOutSamplePhase = atan2f(cimagf(PLLOutSample),crealf(PLLOutSample));
         #endif
         //PLLOutSamplePhase = atan2_approximation1d(cimag(PLLOutSample),creal(PLLOutSample));
         
         //running average of the absolute value of the sample phase
         averagePhase = averagePhase * (1.0 - averagePhaseAlpha) + averagePhaseAlpha * fabsf(PLLOutSamplePhase);
         
         //Calculate Error
         #if USE_ATAN_APPROX==1
         sample_phase = arctan2(cimagf(complexDataIn[idx]),crealf(complexDataIn[idx]));
         #else
         sample_phase = atan2f(cimagf(complexDataIn[idx]),crealf(complexDataIn[idx]));
         #endif
         //sample_phase = atan2_approximation1d(cimag(complexDataIn[idx]),creal(complexDataIn[idx]));
      #else
         t_imag = sin(d_phase);
         t_real = cos(d_phase);
         
         //shift the frequency by the carrier 
         PLLOutSample = complexDataIn[idx] * (t_real+I*-t_imag);
          
         //data bits are in the imaginary part
         realDataOut[idx] = cimag(PLLOutSample);
          
         //calculate phase angle for quality estimation 
         #if USE_ATAN_APPROX==1
         PLLOutSamplePhase = arctan2(cimag(PLLOutSample),creal(PLLOutSample));
         #else
         PLLOutSamplePhase = atan2(cimag(PLLOutSample),creal(PLLOutSample));
         #endif
         //PLLOutSamplePhase = atan2_approximation1d(cimag(PLLOutSample),creal(PLLOutSample));
         
         //running average of the absolute value of the sample phase
         averagePhase = averagePhase * (1.0 - averagePhaseAlpha) + averagePhaseAlpha * fabs(PLLOutSamplePhase);
         
         //Calculate Error
         #if USE_ATAN_APPROX==1
         sample_phase = arctan2(cimag(complexDataIn[idx]),creal(complexDataIn[idx]));
         #else
         sample_phase = atan2(cimag(complexDataIn[idx]),creal(complexDataIn[idx]));
         #endif
         //sample_phase = atan2_approximation1d(cimag(complexDataIn[idx]),creal(complexDataIn[idx]));
      #endif
      
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
         d_phase = d_phase-2.0*M_PI;
      
      while(d_phase < -2*M_PI)
         d_phase = d_phase+2.0*M_PI;
   
      //Limit the frequency
      if(d_freq > d_max_freq)
         d_freq = d_max_freq;
      else if(d_freq < d_min_freq)
         d_freq = d_min_freq;
       
      //d_locksig = d_locksig * (1.0 - d_alpha) + d_alpha*(real(dataStreamIn(idx)) * t_real + imag(dataStreamIn(idx)) * t_imag);
      //d_locksig = d_locksig * (1.0 - lockSigAlpha) + lockSigAlpha*(creal(complexDataIn[idx]) * t_real + cimag(complexDataIn[idx]) * t_imag);
      
      #if USE_FLOATS==1
         complexDataInRe = crealf(complexDataIn[idx]);
         complexDataInIm = cimagf(complexDataIn[idx]);
         magnitude_squared = complexDataInRe*complexDataInRe+complexDataInIm*complexDataInIm;
         invsqrt = Q_rsqrt(magnitude_squared);
      #else
         complexDataInRe = creal(complexDataIn[idx]);
         complexDataInIm = cimag(complexDataIn[idx]);
         magnitude_squared = complexDataInRe*complexDataInRe+complexDataInIm*complexDataInIm;
         invsqrt = Q_rsqrt(magnitude_squared); //what is a float function doing here
      #endif
      
      
      
      
      //invsqrt = 1/sqrt(magnitude_squared);
      
      complexDataInRe *= invsqrt;
      complexDataInIm *= invsqrt;
      
      //printf("%f ",sqrt(complexDataInRe*complexDataInRe+complexDataInIm*complexDataInIm));
      //normalizedDataIn = complexDataIn[idx]/cabsf(complexDataIn[idx]);
      //normalizedDataIn = complexDataIn[idx]*Q_rsqrt(powf(crealf(complexDataIn[idx]),2)+powf(cimagf(complexDataIn[idx]),2));
      //normalizedDataIn = complexDataIn[idx]*Q_rsqrt(pow(creal(complexDataIn[idx]),2)+pow(cimag(complexDataIn[idx]),2));
      //normalizedDataIn = complexDataIn[idx];
      
      //d_locksig = d_locksig * (1.0 - lockSigAlpha) + lockSigAlpha*(creal(normalizedDataIn) * t_real + cimag(normalizedDataIn) * t_imag);
      d_locksig = d_locksig * (1.0 - lockSigAlpha) + lockSigAlpha*(complexDataInRe * t_real + complexDataInIm * t_imag);
      
      if(lockSignalStreamOut != NULL)  
          lockSignalStreamOut[idx] = d_locksig;
      
       //sweep the frequency
       //if(idx > 32000 && d_locksig < d_lock_threshold/15.0 && firstLock == 0)
       //if(d_locksig < d_lock_threshold/15.0 && firstLock == -1)
       //if(d_locksig < d_lock_threshold/3.0 && firstLock == -1)
       //if(fabs(CONST_CNTR_ANGLE-averagePhase) < 0.05 && firstLock == -1)
       
       #if USE_FLOATS==1
          if(fabsf(CONST_CNTR_ANGLE-averagePhase) < 0.05  && firstLock == -1)
             {          
             d_freq  = d_freq + sweep;
             if(d_freq >= d_max_freq)
                sweep = sweep*-1.0;
             else if(d_freq <= d_min_freq)
                sweep = sweep*-1.0;
             else
                {
                if(d_freq >= 0)
                   sweep = fabsf(sweep);
                else
                   sweep = fabsf(sweep) * -1.0;
                }
             }
       #else
          if(fabs(CONST_CNTR_ANGLE-averagePhase) < 0.05  && firstLock == -1)
             {          
             d_freq  = d_freq + sweep;
             if(d_freq >= d_max_freq)
                sweep = sweep*-1.0;
             else if(d_freq <= d_min_freq)
                sweep = sweep*-1.0;
             else
                {
                if(d_freq >= 0)
                   sweep = fabs(sweep);
                else
                   sweep = fabs(sweep) * -1.0;
                }
             }
       #endif
       
      //d_locksig > d_lock_threshold (0.01)
      if(d_locksig > d_lock_threshold && firstLock == -1)
      //if((1.5708-averagePhase) > d_lock_threshold && firstLock == -1)
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